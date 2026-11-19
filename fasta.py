"""Random access to a reference genome FASTA, standard library only.

Plain language
--------------
A reference FASTA is a plain text file that holds the DNA sequence of every
chromosome, written as lines of letters under a ``>name`` title line. It is
large (the human genome is about 3 GB), so we never read it into memory. We
instead read a small side-file - the **index**, suffix ``.fai`` - that records,
for every chromosome, *where in the file its sequence starts* and *how many
letters are on each line*. With those four numbers we can compute the exact
byte position of any base and seek straight to it.

The ``.fai`` format is the samtools faidx format, five tab-separated columns::

    name    length    offset    linebases    linewidth

``offset`` is the byte position of the first base of that chromosome,
``linebases`` the number of bases on a full line, and ``linewidth`` the number
of *bytes* on a full line (``linebases`` plus the newline, which is 1 byte on
Linux/macOS and 2 bytes on Windows). GTT writes an index that samtools can
read, and reads an index that samtools wrote.

Requirement: every sequence line of a record must have the same length, except
the final line of that record. That is what makes the arithmetic possible, and
it is also what ``samtools faidx`` requires.
"""
from __future__ import annotations

import os
from dataclasses import dataclass
from typing import Dict, List, Optional

from .errors import FastaError

_GZIP_MAGIC = b"\x1f\x8b"

#: Bases fetched per read when walking left during left-alignment.
DEFAULT_WINDOW = 64


def _canonical(name: str) -> str:
    """Canonical contig name: drop a ``chr`` prefix, fold ``M`` onto ``MT``."""
    n = name[3:] if name.lower().startswith("chr") else name
    return "MT" if n.upper() == "M" else n


@dataclass(slots=True)
class FaiRecord:
    """One line of a ``.fai`` index."""

    name: str
    length: int
    offset: int
    linebases: int
    linewidth: int

    def byte_offset(self, pos0: int) -> int:
        """Byte offset of the 0-based base ``pos0`` of this contig."""
        return self.offset + (pos0 // self.linebases) * self.linewidth + (pos0 % self.linebases)


def build_fai(fasta_path: str) -> List[FaiRecord]:
    """Scan ``fasta_path`` once and compute its index records.

    Raises :class:`~gtt.errors.FastaError` if the file is gzip-compressed, has
    no records, has an empty record, or has ragged line lengths inside a
    record (the one layout the ``.fai`` arithmetic cannot describe).
    """
    if not os.path.exists(fasta_path):
        raise FastaError(f"reference FASTA does not exist: {fasta_path}")
    with open(fasta_path, "rb") as fh:
        if fh.read(2) == _GZIP_MAGIC:
            raise FastaError(
                f"{fasta_path} is gzip-compressed. Indexed random access needs an "
                "uncompressed FASTA (or a bgzip-compressed one with a .gzi index, which "
                "GTT does not read). Decompress it first, for example with "
                "`bgzip -d reference.fa.gz` or 7-Zip."
            )

    records: List[FaiRecord] = []
    name: Optional[str] = None
    offset = 0
    length = 0
    linebases = 0
    linewidth = 0
    short_line_seen = False
    line_no = 0

    def _close_record():
        if name is None:
            return
        if length == 0:
            raise FastaError(f"{fasta_path}: contig {name!r} has no sequence")
        records.append(FaiRecord(name, length, offset, linebases, linewidth))

    with open(fasta_path, "rb") as fh:
        byte_pos = 0
        for raw in fh:
            line_no += 1
            start_of_line = byte_pos
            byte_pos += len(raw)
            stripped = raw.rstrip(b"\r\n")
            if raw.startswith(b">"):
                _close_record()
                name = stripped[1:].split()[0].decode("ascii", "replace") if len(stripped) > 1 else ""
                if not name:
                    raise FastaError(f"{fasta_path}:{line_no}: '>' line has no sequence name")
                offset = byte_pos
                length = 0
                linebases = 0
                linewidth = 0
                short_line_seen = False
                continue
            if name is None:
                if not stripped:
                    continue
                raise FastaError(f"{fasta_path}:{line_no}: sequence data before any '>' line")
            if not stripped:
                # A blank line inside a record ends it for indexing purposes.
                short_line_seen = True
                continue
            if short_line_seen:
                raise FastaError(
                    f"{fasta_path}:{line_no}: contig {name!r} has a short line before its last "
                    "line. Every sequence line of a record must be the same length (except the "
                    "last) for an .fai index to be possible. Re-wrap the FASTA, for example with "
                    "`samtools faidx` after `seqkit seq -w 60`."
                )
            if linebases == 0:
                linebases = len(stripped)
                linewidth = len(raw)
            elif len(stripped) < linebases:
                short_line_seen = True
            elif len(stripped) > linebases:
                raise FastaError(
                    f"{fasta_path}:{line_no}: line is longer ({len(stripped)}) than the first "
                    f"line of contig {name!r} ({linebases}); ragged FASTA cannot be indexed"
                )
            length += len(stripped)
            _ = start_of_line
        _close_record()

    if not records:
        raise FastaError(f"{fasta_path}: no FASTA records found (no '>' line)")
    seen = set()
    for rec in records:
        if rec.name in seen:
            raise FastaError(f"{fasta_path}: contig {rec.name!r} appears more than once")
        seen.add(rec.name)
    return records


def write_fai(records: List[FaiRecord], fai_path: str) -> None:
    """Write ``records`` to ``fai_path`` in samtools faidx format."""
    with open(fai_path, "w", newline="\n", encoding="ascii") as fh:
        for r in records:
            fh.write(f"{r.name}\t{r.length}\t{r.offset}\t{r.linebases}\t{r.linewidth}\n")


def read_fai(fai_path: str) -> List[FaiRecord]:
    """Parse an existing ``.fai`` index."""
    records: List[FaiRecord] = []
    with open(fai_path, "rt", encoding="ascii", errors="replace") as fh:
        for i, line in enumerate(fh, start=1):
            line = line.rstrip("\r\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 5:
                raise FastaError(f"{fai_path}:{i}: expected 5 tab-separated columns, got {len(parts)}")
            try:
                records.append(FaiRecord(parts[0], int(parts[1]), int(parts[2]),
                                         int(parts[3]), int(parts[4])))
            except ValueError:
                raise FastaError(f"{fai_path}:{i}: non-integer value in index line") from None
    if not records:
        raise FastaError(f"{fai_path}: index file is empty")
    return records


class IndexedFasta:
    """Read-only random access to a reference FASTA.

    Parameters
    ----------
    path:
        The uncompressed reference FASTA.
    fai_path:
        Index location. Defaults to ``path + '.fai'``.
    create_index:
        Build the index if it is missing. If the directory cannot be written
        to, the index is kept in memory instead and :attr:`index_written` is
        False.

    Contig names are matched leniently: ``chr20``, ``20``, ``chrM`` and ``MT``
    all find each other, because VCF files and reference files disagree about
    this constantly and the disagreement is never biological.
    """

    def __init__(self, path: str, fai_path: Optional[str] = None, create_index: bool = True):
        self.path = path
        self.fai_path = fai_path or (path + ".fai")
        self.index_written = False
        self.index_created = False
        self.notes: List[str] = []
        if os.path.exists(self.fai_path):
            records = read_fai(self.fai_path)
            stat_fa = os.path.getmtime(path) if os.path.exists(path) else 0
            if os.path.getmtime(self.fai_path) < stat_fa:
                self.notes.append(
                    f"{self.fai_path} is older than {os.path.basename(path)}; it is used as "
                    "found. Delete it to force a rebuild if the FASTA changed."
                )
        else:
            records = build_fai(path)
            self.index_created = True
            if create_index:
                try:
                    write_fai(records, self.fai_path)
                    self.index_written = True
                    self.notes.append(f"created index {self.fai_path}")
                except OSError as exc:
                    self.notes.append(
                        f"could not write {self.fai_path} ({exc}); index held in memory only"
                    )
        self.records: Dict[str, FaiRecord] = {r.name: r for r in records}
        self._alias: Dict[str, str] = {}
        for r in records:
            self._alias.setdefault(_canonical(r.name), r.name)
        self._fh = open(path, "rb")
        self.n_fetches = 0
        self.bases_fetched = 0

    # -- contig names -----------------------------------------------------
    @property
    def contigs(self) -> Dict[str, int]:
        """Mapping of contig name to length, in file order."""
        return {name: rec.length for name, rec in self.records.items()}

    def resolve(self, contig: str) -> Optional[str]:
        """Return the FASTA's own name for ``contig``, or None if absent."""
        if contig in self.records:
            return contig
        return self._alias.get(_canonical(contig))

    def length(self, contig: str) -> Optional[int]:
        name = self.resolve(contig)
        return self.records[name].length if name else None

    # -- sequence ---------------------------------------------------------
    def fetch(self, contig: str, start: int, end: int) -> str:
        """Return bases ``[start, end)`` of ``contig``, 0-based, uppercased.

        >>> # bases 0..3 of chr1 are the first three letters of that sequence
        >>> # fa.fetch("chr1", 0, 3)
        """
        name = self.resolve(contig)
        if name is None:
            raise FastaError(f"contig {contig!r} is not in {os.path.basename(self.path)}")
        rec = self.records[name]
        if start < 0 or end > rec.length or end < start:
            raise FastaError(
                f"requested {contig}:{start}-{end} (0-based) outside {name} "
                f"which is {rec.length} bases long"
            )
        if end == start:
            return ""
        b0 = rec.byte_offset(start)
        b1 = rec.byte_offset(end - 1) + 1
        self._fh.seek(b0)
        raw = self._fh.read(b1 - b0)
        seq = raw.replace(b"\r", b"").replace(b"\n", b"").decode("ascii", "replace").upper()
        if len(seq) != end - start:
            raise FastaError(
                f"{self.path}: index does not match the file at {name}:{start}-{end} "
                f"(expected {end - start} bases, read {len(seq)}). The .fai index is stale; "
                "delete it and let GTT rebuild it."
            )
        self.n_fetches += 1
        self.bases_fetched += len(seq)
        return seq

    def base(self, contig: str, pos1: int) -> str:
        """Return the single base at 1-based position ``pos1``."""
        return self.fetch(contig, pos1 - 1, pos1)

    def close(self):
        try:
            self._fh.close()
        except Exception:
            pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def summary(self) -> Dict[str, object]:
        return {
            "path": os.path.basename(self.path),
            "fai": os.path.basename(self.fai_path),
            "index_created_by_gtt": self.index_created,
            "index_written_to_disk": self.index_written,
            "n_contigs": len(self.records),
            "notes": list(self.notes),
        }
