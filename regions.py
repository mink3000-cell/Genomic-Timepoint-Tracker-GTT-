"""BED region files: parsing, merging, intersection and membership tests.

Plain language
--------------
A **BED file** lists stretches of a chromosome, one per line::

    chr20    100    200

That line means "chromosome 20, from base 101 to base 200 inclusive". BED
counts from 0 and leaves the end out (*half-open*), while VCF counts from 1 and
includes it. The conversion is therefore ``bed_start = vcf_pos - 1``. Getting
this wrong by one base is the classic bioinformatics bug, so it is done in
exactly one place here: :meth:`RegionSet.contains`.

Why GTT needs this
------------------
A variant caller only makes a statement about the parts of the genome it could
actually examine. Outside those parts the file says nothing - and *nothing* is
not the same as *no variant*. When two files have different confident regions,
every position confident in one and not in the other produces a difference that
is purely about coverage. Restricting the comparison to the **intersection** of
both region sets removes that artefact. Positions outside the intersection are
reported by GTT as **not interrogated**, never as absent.
"""
from __future__ import annotations

import bisect
import gzip
import io
import os
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

from .errors import BedFormatError

_GZIP_MAGIC = b"\x1f\x8b"


def _canonical(name: str) -> str:
    """Canonical contig name: drop a ``chr`` prefix, fold ``M`` onto ``MT``."""
    n = name[3:] if name.lower().startswith("chr") else name
    return "MT" if n.upper() == "M" else n


def _open_text(path: str):
    if not os.path.exists(path):
        raise BedFormatError(f"BED file does not exist: {path}")
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == _GZIP_MAGIC:
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


def merge_intervals(intervals: Iterable[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Sort and merge overlapping or touching half-open intervals."""
    items = sorted(i for i in intervals if i[1] > i[0])
    out: List[Tuple[int, int]] = []
    for start, end in items:
        if out and start <= out[-1][1]:
            if end > out[-1][1]:
                out[-1] = (out[-1][0], end)
        else:
            out.append((start, end))
    return out


def intersect_intervals(a: Sequence[Tuple[int, int]],
                        b: Sequence[Tuple[int, int]]) -> List[Tuple[int, int]]:
    """Intersect two already-merged, sorted interval lists."""
    out: List[Tuple[int, int]] = []
    i = j = 0
    while i < len(a) and j < len(b):
        start = max(a[i][0], b[j][0])
        end = min(a[i][1], b[j][1])
        if end > start:
            out.append((start, end))
        if a[i][1] <= b[j][1]:
            i += 1
        else:
            j += 1
    return out


class RegionSet:
    """A set of genomic intervals, held per contig as merged half-open spans.

    Contig names are stored canonically (``chr20`` and ``20`` are the same
    contig), because a BED file and a VCF file disagreeing about the ``chr``
    prefix is a naming accident, not a biological statement.
    """

    def __init__(self, intervals: Optional[Dict[str, List[Tuple[int, int]]]] = None,
                 label: str = "regions", sources: Optional[List[str]] = None):
        self.label = label
        self.sources: List[str] = list(sources or [])
        self.by_contig: Dict[str, List[Tuple[int, int]]] = {}
        self.names_seen: Dict[str, str] = {}
        for contig, spans in (intervals or {}).items():
            canon = _canonical(contig)
            self.names_seen.setdefault(canon, contig)
            self.by_contig[canon] = merge_intervals(
                list(self.by_contig.get(canon, [])) + list(spans))
        self._starts: Dict[str, List[int]] = {
            c: [s for s, _ in spans] for c, spans in self.by_contig.items()}

    # -- construction -----------------------------------------------------
    @classmethod
    def from_bed(cls, path: str, label: Optional[str] = None) -> "RegionSet":
        """Parse a BED file (plain or gzip). Comment/track/browser lines are skipped."""
        intervals: Dict[str, List[Tuple[int, int]]] = {}
        names_seen: Dict[str, str] = {}
        n_lines = 0
        stream = _open_text(path)
        try:
            for line_no, line in enumerate(stream, start=1):
                line = line.rstrip("\r\n")
                if not line.strip():
                    continue
                if line.startswith(("#", "track", "browser")):
                    continue
                parts = line.split("\t") if "\t" in line else line.split()
                if len(parts) < 3:
                    raise BedFormatError(
                        f"{path}:{line_no}: BED needs at least 3 columns "
                        f"(chrom, start, end), found {len(parts)}"
                    )
                contig = parts[0]
                try:
                    start, end = int(parts[1]), int(parts[2])
                except ValueError:
                    raise BedFormatError(
                        f"{path}:{line_no}: start/end are not integers: {parts[1]!r}, {parts[2]!r}"
                    ) from None
                if start < 0:
                    raise BedFormatError(f"{path}:{line_no}: negative start {start}")
                if end < start:
                    raise BedFormatError(
                        f"{path}:{line_no}: end {end} is before start {start}")
                n_lines += 1
                canon = _canonical(contig)
                names_seen.setdefault(canon, contig)
                intervals.setdefault(canon, []).append((start, end))
        finally:
            stream.close()
        if not intervals:
            raise BedFormatError(f"{path}: no usable BED records found")
        rs = cls(intervals, label=label or os.path.basename(path), sources=[path])
        rs.names_seen.update(names_seen)
        rs.n_input_lines = n_lines
        return rs

    def intersect(self, other: "RegionSet", label: Optional[str] = None) -> "RegionSet":
        """Return the intersection of two region sets (contig-wise)."""
        out: Dict[str, List[Tuple[int, int]]] = {}
        for contig, spans in self.by_contig.items():
            if contig in other.by_contig:
                shared = intersect_intervals(spans, other.by_contig[contig])
                if shared:
                    out[contig] = shared
        rs = RegionSet(out, label=label or f"{self.label} AND {other.label}",
                       sources=self.sources + other.sources)
        rs.names_seen.update(self.names_seen)
        return rs

    # -- queries ----------------------------------------------------------
    def contains(self, contig: str, pos1: int) -> bool:
        """Is the 1-based position ``pos1`` inside this region set?

        The only place where the BED/VCF coordinate conversion happens.
        """
        canon = _canonical(contig)
        spans = self.by_contig.get(canon)
        if not spans:
            return False
        p = pos1 - 1
        i = bisect.bisect_right(self._starts[canon], p) - 1
        return i >= 0 and spans[i][0] <= p < spans[i][1]

    def contains_span(self, contig: str, pos1: int, ref_len: int) -> bool:
        """Is every base of ``[pos1, pos1+ref_len)`` inside this region set?"""
        canon = _canonical(contig)
        spans = self.by_contig.get(canon)
        if not spans:
            return False
        p = pos1 - 1
        i = bisect.bisect_right(self._starts[canon], p) - 1
        return i >= 0 and spans[i][0] <= p and p + max(ref_len, 1) <= spans[i][1]

    def has_contig(self, contig: str) -> bool:
        return _canonical(contig) in self.by_contig

    @property
    def contigs(self) -> List[str]:
        return list(self.by_contig)

    @property
    def n_intervals(self) -> int:
        return sum(len(v) for v in self.by_contig.values())

    @property
    def total_bases(self) -> int:
        return sum(end - start for spans in self.by_contig.values() for start, end in spans)

    def bases_on(self, contig: str) -> int:
        spans = self.by_contig.get(_canonical(contig), [])
        return sum(end - start for start, end in spans)

    def summary(self) -> Dict[str, object]:
        return {
            "label": self.label,
            "sources": [os.path.basename(s) for s in self.sources],
            "n_contigs": len(self.by_contig),
            "n_intervals": self.n_intervals,
            "total_bases": self.total_bases,
        }

    def __repr__(self):
        return (f"RegionSet({self.label!r}, {len(self.by_contig)} contigs, "
                f"{self.n_intervals} intervals, {self.total_bases} bases)")


def build_comparison_regions(regions: Optional[str] = None,
                             regions_baseline: Optional[str] = None,
                             regions_followup: Optional[str] = None
                             ) -> Tuple[Optional[RegionSet], Dict[str, object]]:
    """Combine up to three BED files into one comparison region set.

    ``regions`` applies to both inputs; ``regions_baseline`` and
    ``regions_followup`` describe one input each. When more than one is given
    the comparison region is their **intersection**: a position must be
    interrogated in every supplied region file to be compared.

    Returns ``(region_set_or_None, provenance_dict)``.
    """
    parts: List[RegionSet] = []
    prov: Dict[str, object] = {}
    for key, path in (("regions", regions), ("regions_baseline", regions_baseline),
                      ("regions_followup", regions_followup)):
        if path:
            rs = RegionSet.from_bed(path, label=key)
            parts.append(rs)
            prov[key] = rs.summary()
    if not parts:
        return None, prov
    combined = parts[0]
    for rs in parts[1:]:
        combined = combined.intersect(rs)
    combined.label = "comparison regions"
    prov["comparison_regions"] = combined.summary()
    prov["combination"] = (
        "intersection of " + ", ".join(os.path.basename(s) for s in combined.sources)
        if len(parts) > 1 else os.path.basename(combined.sources[0])
    )
    return combined, prov
