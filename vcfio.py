"""Streaming VCF input: header metadata and record decomposition.

The reader never holds the whole file in memory; it yields one
:class:`~gtt.model.AlleleCall` per (record, ALT allele) pair.
"""
from __future__ import annotations

import gzip
import io
import os
import re
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Tuple

from .errors import VcfFormatError
from .model import (
    AlleleCall,
    GT_HEMI,
    GT_HET,
    GT_HOM_ALT,
    GT_HOM_REF,
    GT_NO_CALL,
    GT_NO_GT_FIELD,
    GT_NO_SAMPLE,
    GT_OTHER_ALT,
)
from .normalize import (
    LA_AT_CONTIG_START,
    LA_CONTIG_ABSENT,
    LA_REF_MISMATCH,
    LA_SHIFT_LIMIT,
    LA_SHIFTED,
    LA_UNCHANGED,
    MAX_LEFT_SHIFT,
    NON_VARIANT_ALTS,
    classify_variant,
    is_plain_dna,
    is_symbolic,
    left_align_allele,
    normalize_allele,
)

MANDATORY_COLUMNS = ["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]
_META_KV_RE = re.compile(r'(\w+)=(?:"((?:[^"\\]|\\.)*)"|([^,>]*))')
_GZIP_MAGIC = b"\x1f\x8b"


def _parse_meta_fields(text: str) -> Dict[str, str]:
    """Parse the ``<ID=x,Number=1,...>`` payload of a structured meta line."""
    inner = text.strip()
    if inner.startswith("<") and inner.endswith(">"):
        inner = inner[1:-1]
    out = {}
    for m in _META_KV_RE.finditer(inner):
        out[m.group(1)] = m.group(2) if m.group(2) is not None else m.group(3)
    return out


def open_vcf(path: str):
    """Open a plain or (b)gzip-compressed VCF as a text stream, by magic bytes."""
    if not os.path.exists(path):
        raise VcfFormatError("file does not exist", path=path)
    if os.path.getsize(path) == 0:
        raise VcfFormatError("file is empty", path=path)
    with open(path, "rb") as fh:
        magic = fh.read(2)
    if magic == _GZIP_MAGIC:
        return io.TextIOWrapper(gzip.open(path, "rb"), encoding="utf-8", errors="replace")
    return open(path, "rt", encoding="utf-8", errors="replace")


@dataclass
class VcfHeader:
    """Structured view of a VCF header."""

    path: str
    fileformat: Optional[str] = None
    reference: Optional[str] = None
    contigs: Dict[str, Optional[int]] = field(default_factory=dict)
    contig_assembly: Dict[str, str] = field(default_factory=dict)
    filters: Dict[str, str] = field(default_factory=dict)
    info: Dict[str, Dict[str, str]] = field(default_factory=dict)
    formats: Dict[str, Dict[str, str]] = field(default_factory=dict)
    samples: List[str] = field(default_factory=list)
    other_meta: List[str] = field(default_factory=list)
    n_meta_lines: int = 0
    synthetic: bool = False
    synthetic_note: Optional[str] = None

    @property
    def chr_prefixed(self) -> Optional[bool]:
        """True if contigs are ``chr``-prefixed, False if not, None if unknown."""
        if not self.contigs:
            return None
        names = list(self.contigs)
        return all(n.startswith("chr") for n in names)

    @property
    def n_samples(self) -> int:
        return len(self.samples)

    def summary(self) -> Dict[str, object]:
        return {
            "path": os.path.basename(self.path),
            "fileformat": self.fileformat,
            "reference": self.reference,
            "n_contigs": len(self.contigs),
            "samples": list(self.samples),
            "n_filters": len(self.filters),
            "chr_prefixed": self.chr_prefixed,
            "synthetic": self.synthetic,
        }


def read_header(path: str) -> VcfHeader:
    """Read only the header of ``path``."""
    stream = open_vcf(path)
    try:
        hdr, _, _ = _consume_header(stream, path)
    finally:
        stream.close()
    return hdr


def _consume_header(stream, path) -> Tuple[VcfHeader, List[str], int]:
    """Read meta + column lines; return (header, column_names, lines_consumed)."""
    hdr = VcfHeader(path=path)
    columns: List[str] = []
    line_no = 0
    for line in stream:
        line_no += 1
        line = line.rstrip("\n").rstrip("\r")
        if line.startswith("##"):
            hdr.n_meta_lines += 1
            key, _, value = line[2:].partition("=")
            if key == "fileformat":
                hdr.fileformat = value
            elif key == "reference":
                hdr.reference = value
            elif key == "contig":
                f = _parse_meta_fields(value)
                cid = f.get("ID")
                if cid:
                    length = f.get("length")
                    hdr.contigs[cid] = int(length) if length and length.isdigit() else None
                    if f.get("assembly"):
                        hdr.contig_assembly[cid] = f["assembly"]
            elif key == "FILTER":
                f = _parse_meta_fields(value)
                if f.get("ID"):
                    hdr.filters[f["ID"]] = f.get("Description", "")
            elif key == "INFO":
                f = _parse_meta_fields(value)
                if f.get("ID"):
                    hdr.info[f["ID"]] = f
            elif key == "FORMAT":
                f = _parse_meta_fields(value)
                if f.get("ID"):
                    hdr.formats[f["ID"]] = f
            elif key == "gtt_synthetic":
                hdr.synthetic = True
                hdr.synthetic_note = value
            else:
                hdr.other_meta.append(line)
        elif line.startswith("#CHROM"):
            columns = line[1:].split("\t")
            break
        elif not line.strip():
            continue
        else:
            raise VcfFormatError(
                "data line encountered before the #CHROM header line", path=path, line_no=line_no
            )
    if not columns:
        raise VcfFormatError("no #CHROM header line found", path=path)
    if columns[: len(MANDATORY_COLUMNS)] != MANDATORY_COLUMNS:
        raise VcfFormatError(
            f"#CHROM line does not begin with the 8 mandatory VCF columns, found: {columns[:8]}",
            path=path,
            line_no=line_no,
        )
    if len(columns) == 9:
        raise VcfFormatError(
            "#CHROM line has a FORMAT column but no sample columns", path=path, line_no=line_no
        )
    if len(columns) > 9:
        hdr.samples = columns[9:]
    if hdr.fileformat is None:
        raise VcfFormatError("missing mandatory ##fileformat meta line", path=path)
    return hdr, columns, line_no


def _parse_gt(gt_raw: str, alt_index: int) -> Tuple[str, Optional[int], bool]:
    """Project a GT string onto one ALT index.

    Returns ``(gt_status, allele_copies, phased)``.
    """
    phased = "|" in gt_raw
    tokens = re.split(r"[/|]", gt_raw)
    if not tokens or all(t == "." or t == "" for t in tokens):
        return GT_NO_CALL, None, phased
    called = []
    for t in tokens:
        if t == "." or t == "":
            continue
        if not t.isdigit():
            return GT_NO_CALL, None, phased
        called.append(int(t))
    if not called:
        return GT_NO_CALL, None, phased
    copies = called.count(alt_index)
    if copies == 0:
        if all(a == 0 for a in called):
            return GT_HOM_REF, 0, phased
        return GT_OTHER_ALT, 0, phased
    if len(called) == 1:
        return GT_HEMI, copies, phased
    if copies == len(called):
        return GT_HOM_ALT, copies, phased
    return GT_HET, copies, phased


@dataclass
class ReadStats:
    """Counters describing one pass over a VCF."""

    records: int = 0
    alleles_emitted: int = 0
    multiallelic_records: int = 0
    skipped_star: int = 0
    skipped_no_alt: int = 0
    skipped_symbolic: int = 0
    skipped_not_pass: int = 0
    symbolic_emitted: int = 0
    not_pass_records: int = 0
    no_call_records: int = 0
    hom_ref_records: int = 0
    # -- v0.2.0: reference-based normalization --
    left_aligned_unchanged: int = 0
    left_aligned_shifted: int = 0
    ref_mismatches: int = 0
    alleles_on_contigs_not_in_reference: int = 0
    truncated_at_contig_start: int = 0
    shift_limit_reached: int = 0
    # -- v0.2.0: region restriction --
    not_interrogated: int = 0
    contigs_seen: List[str] = field(default_factory=list)
    contigs_not_in_reference: List[str] = field(default_factory=list)
    ref_mismatch_examples: List[Dict[str, object]] = field(default_factory=list)
    not_interrogated_by_type: Dict[str, int] = field(default_factory=dict)
    #: Which v0.2.0 features were active. When both are off the counter set -
    #: and therefore the JSON summary - is exactly the v0.1.0 one.
    used_reference: bool = False
    used_regions: bool = False

    #: Counters reported only when ``--reference`` was given.
    REFERENCE_KEYS = ("left_aligned_unchanged", "left_aligned_shifted", "ref_mismatches",
                      "alleles_on_contigs_not_in_reference", "truncated_at_contig_start",
                      "shift_limit_reached")
    #: Counters reported only when ``--regions`` was given.
    REGION_KEYS = ("not_interrogated",)

    def as_dict(self) -> Dict[str, object]:
        skip = {"contigs_seen", "contigs_not_in_reference", "ref_mismatch_examples",
                "not_interrogated_by_type", "used_reference", "used_regions"}
        if not self.used_reference:
            skip.update(self.REFERENCE_KEYS)
        if not self.used_regions:
            skip.update(self.REGION_KEYS)
        d = {k: getattr(self, k) for k in self.__dataclass_fields__ if k not in skip}
        d["n_contigs_seen"] = len(self.contigs_seen)
        if self.used_reference:
            d["contigs_not_in_reference"] = list(self.contigs_not_in_reference)
        if self.used_regions:
            d["not_interrogated_by_type"] = dict(self.not_interrogated_by_type)
        return d


class VcfReader:
    """Streaming reader yielding one :class:`AlleleCall` per record/ALT pair.

    Parameters
    ----------
    path:
        Path to a plain or (b)gzip-compressed VCF.
    sample:
        Sample name to take genotypes from. Defaults to the first sample column.
    skip_star:
        Drop ``*`` ALT alleles (upstream-deletion placeholders, not independent
        variants). Counted in :attr:`stats.skipped_star`.
    skip_symbolic:
        Drop symbolic/breakend ALTs, which trimming cannot normalize.
    pass_only:
        Treat only records with FILTER ``PASS``/``.`` as reported.
    fasta:
        An :class:`~gtt.fasta.IndexedFasta`. When given, alleles are fully
        left-aligned after trimming and every REF is checked against the
        reference sequence. When None, behaviour is exactly v0.1.0's.
    regions:
        A :class:`~gtt.regions.RegionSet`. Alleles outside it are still yielded
        but carry ``interrogated=False``, so the caller can report them as
        *not interrogated* rather than as absent.
    region_mode:
        ``start`` (default) tests the normalized POS only - the convention used
        by ``bcftools -T``. ``span`` requires every base of REF to be inside a
        single interval, which is stricter at region edges.
    """

    def __init__(self, path, sample=None, skip_star=True, skip_symbolic=False, pass_only=False,
                 fasta=None, regions=None, region_mode="start",
                 max_left_shift=MAX_LEFT_SHIFT, max_ref_mismatch_examples=10000):
        self.path = path
        self.requested_sample = sample
        self.skip_star = skip_star
        self.skip_symbolic = skip_symbolic
        self.pass_only = pass_only
        self.fasta = fasta
        self.regions = regions
        if region_mode not in ("start", "span"):
            raise ValueError(f"region_mode must be 'start' or 'span', got {region_mode!r}")
        self.region_mode = region_mode
        self.max_left_shift = max_left_shift
        self.max_ref_mismatch_examples = max_ref_mismatch_examples
        self.stats = ReadStats(used_reference=fasta is not None, used_regions=regions is not None)
        self._stream = open_vcf(path)
        try:
            self.header, self._columns, self._line_no = _consume_header(self._stream, path)
        except Exception:
            self.close()
            raise
        if sample is not None:
            if sample not in self.header.samples:
                self.close()
                raise VcfFormatError(
                    f"sample {sample!r} not found; file has {self.header.samples}", path=path
                )
            self.sample = sample
        else:
            self.sample = self.header.samples[0] if self.header.samples else None
        self._sample_col = 9 + self.header.samples.index(self.sample) if self.sample else None
        self._n_cols = len(self._columns)

    def close(self):
        try:
            self._stream.close()
        except Exception:
            pass

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()

    def __iter__(self) -> Iterator[AlleleCall]:
        try:
            yield from self._iter_records()
        finally:
            self.close()

    def _count_left_align(self, contig, orig_pos, orig_ref, orig_alt,
                          npos, nref, nalt, status, ref_seen):
        """Tally one left-alignment outcome; keep REF mismatches visible."""
        st = self.stats
        if status == LA_SHIFTED:
            st.left_aligned_shifted += 1
        elif status == LA_UNCHANGED:
            st.left_aligned_unchanged += 1
        elif status == LA_AT_CONTIG_START:
            st.truncated_at_contig_start += 1
        elif status == LA_SHIFT_LIMIT:
            st.shift_limit_reached += 1
        elif status == LA_CONTIG_ABSENT:
            st.alleles_on_contigs_not_in_reference += 1
            if contig not in st.contigs_not_in_reference:
                st.contigs_not_in_reference.append(contig)
        elif status == LA_REF_MISMATCH:
            st.ref_mismatches += 1
            if len(st.ref_mismatch_examples) < self.max_ref_mismatch_examples:
                st.ref_mismatch_examples.append({
                    "contig": contig, "pos": orig_pos,
                    "ref_in_file": orig_ref, "alt_in_file": orig_alt,
                    "ref_in_fasta": ref_seen if ref_seen is not None
                    else "(position outside contig)",
                    "line_no": self._line_no,
                })

    def _iter_records(self) -> Iterator[AlleleCall]:
        st = self.stats
        seen = set()
        for line in self._stream:
            self._line_no += 1
            if not line or line[0] == "#":
                continue
            line = line.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 8:
                raise VcfFormatError(
                    f"record has {len(fields)} columns, fewer than the 8 mandatory ones",
                    path=self.path, line_no=self._line_no,
                )
            if len(fields) != self._n_cols:
                raise VcfFormatError(
                    f"record has {len(fields)} columns but the header declares {self._n_cols}",
                    path=self.path, line_no=self._line_no,
                )
            contig, pos_raw, var_id, ref, alt_raw, qual, filt = fields[:7]
            try:
                pos = int(pos_raw)
            except ValueError:
                raise VcfFormatError(
                    f"POS {pos_raw!r} is not an integer", path=self.path, line_no=self._line_no
                ) from None
            if pos < 1:
                raise VcfFormatError(
                    f"POS {pos} is not 1-based positive", path=self.path, line_no=self._line_no
                )
            ref = ref.upper()
            if not ref or not is_plain_dna(ref):
                raise VcfFormatError(
                    f"REF {ref!r} is not a plain [ACGTN] sequence",
                    path=self.path, line_no=self._line_no,
                )
            st.records += 1
            if contig not in seen:
                seen.add(contig)
                st.contigs_seen.append(contig)

            filters = tuple(f for f in filt.split(";") if f) if filt not in (".", "") else ()
            is_pass = filt in ("PASS", ".", "")
            if not is_pass:
                st.not_pass_records += 1
            if self.pass_only and not is_pass:
                st.skipped_not_pass += 1
                continue

            alts = alt_raw.split(",")
            n_alts = len([a for a in alts if a not in (".",)])
            if len(alts) > 1:
                st.multiallelic_records += 1

            gt_raw = None
            if self._sample_col is not None and len(fields) > self._sample_col:
                fmt_keys = fields[8].split(":") if len(fields) > 8 else []
                sample_vals = fields[self._sample_col].split(":")
                if "GT" in fmt_keys:
                    idx = fmt_keys.index("GT")
                    gt_raw = sample_vals[idx] if idx < len(sample_vals) else "./."

            record_no_call = False
            record_hom_ref = False
            for i, alt in enumerate(alts, start=1):
                alt_u = alt.upper()
                if alt_u == ".":
                    st.skipped_no_alt += 1
                    continue
                if alt_u == "*":
                    st.skipped_star += 1
                    if self.skip_star:
                        continue
                symbolic = is_symbolic(alt_u)
                if symbolic:
                    if self.skip_symbolic:
                        st.skipped_symbolic += 1
                        continue
                    st.symbolic_emitted += 1

                npos, nref, nalt = normalize_allele(pos, ref, alt_u)
                la_status = ""
                if self.fasta is not None:
                    npos, nref, nalt, la_status, ref_seen = left_align_allele(
                        self.fasta, contig, npos, nref, nalt, max_shift=self.max_left_shift)
                    self._count_left_align(contig, pos, ref, alt_u, npos, nref, nalt,
                                           la_status, ref_seen)
                interrogated = True
                if self.regions is not None:
                    if self.region_mode == "span":
                        interrogated = self.regions.contains_span(contig, npos, len(nref))
                    else:
                        interrogated = self.regions.contains(contig, npos)
                    if not interrogated:
                        st.not_interrogated += 1
                        vt = classify_variant(nref, nalt)
                        st.not_interrogated_by_type[vt] = \
                            st.not_interrogated_by_type.get(vt, 0) + 1
                if gt_raw is None:
                    gt_status = GT_NO_SAMPLE if self.sample is None else GT_NO_GT_FIELD
                    copies, phased = None, False
                else:
                    gt_status, copies, phased = _parse_gt(gt_raw, i)
                if gt_status == GT_NO_CALL:
                    record_no_call = True
                elif gt_status == GT_HOM_REF:
                    record_hom_ref = True

                st.alleles_emitted += 1
                yield AlleleCall(
                    contig=contig, pos=npos, ref=nref, alt=nalt,
                    orig_pos=pos, orig_ref=ref, orig_alt=alt_u,
                    var_id=var_id, qual=qual, filters=filters, is_pass=is_pass,
                    gt=gt_raw, gt_status=gt_status, allele_copies=copies, phased=phased,
                    multiallelic=len(alts) > 1, n_alts=n_alts, alt_index=i,
                    var_type=classify_variant(nref, nalt), symbolic=symbolic,
                    line_no=self._line_no, sample=self.sample,
                    interrogated=interrogated, la_status=la_status,
                )
            if record_no_call:
                st.no_call_records += 1
            if record_hom_ref:
                st.hom_ref_records += 1
