"""Baseline-versus-follow-up comparison engine.

The engine compares *normalized allele keys* ``(contig, pos, ref, alt)``.
Genotype is not part of the key: a variant reported het in one file and hom in
the other is a shared variant with a reported genotype difference, which is a
different statement from presence or absence.

Memory model: one contig block from each file is held at a time, so peak memory
scales with the largest chromosome, not the genome. Records must be grouped by
contig (the two files may list contigs in different orders); a contig that
reappears after another contig has begun raises
:class:`~gtt.errors.UnsortedInputError`.
"""
from __future__ import annotations

import csv
import os
from dataclasses import dataclass, field
from typing import Dict, Iterator, List, Optional, Tuple

from .errors import UnsortedInputError
from .fasta import IndexedFasta
from .model import AlleleCall, GT_HOM_REF, GT_NO_CALL, format_key
from .normalize import LEVEL_LEFT_ALIGNED, LEVEL_TRIM_ONLY, NORMALIZATION_LEVEL
from .regions import build_comparison_regions
from .validate import CompatibilityReport, check_compatibility, require_compatible
from .vcfio import VcfReader

BASELINE = "baseline"
FOLLOWUP = "followup"

# Why the other file does not report an allele. Never implies biological change.
CTX_NO_RECORD = "no_record_at_site"
CTX_OTHER_ALLELE = "site_record_other_allele"
CTX_FILTERED = "allele_present_but_filtered"
CTX_NO_CALL = "allele_present_but_no_call"
CTX_HOM_REF = "allele_present_but_hom_ref"

CONTEXT_WORDING = {
    CTX_NO_RECORD: "no record at this position in the {other} file",
    CTX_OTHER_ALLELE: "a record exists at this position in the {other} file but reports a different allele",
    CTX_FILTERED: "this allele is present in the {other} file but its record is filtered ({detail})",
    CTX_NO_CALL: "this allele is present in the {other} file but the genotype is a no-call",
    CTX_HOM_REF: "this allele is present in the {other} file but the genotype is homozygous reference",
}

DIFF_COLUMNS = [
    "variant", "contig", "pos", "ref", "alt", "var_type", "reported_in",
    "as_written", "id", "filter", "genotype", "gt_status", "multiallelic",
    "other_file_context", "context_detail", "line_no", "statement",
]
SHARED_COLUMNS = [
    "variant", "contig", "pos", "ref", "alt", "var_type",
    "baseline_as_written", "followup_as_written",
    "baseline_filter", "followup_filter", "baseline_genotype", "followup_genotype",
    "baseline_gt_status", "followup_gt_status",
    "baseline_allele_copies", "followup_allele_copies", "genotype_agrees",
    "baseline_line_no", "followup_line_no",
]
#: Alleles outside the interrogated regions. Listed, never counted as absent.
NOT_INTERROGATED_COLUMNS = [
    "variant", "contig", "pos", "ref", "alt", "var_type", "reported_in",
    "as_written", "filter", "genotype", "line_no", "statement",
]
#: Positions where the file's REF disagrees with the reference FASTA.
REF_MISMATCH_COLUMNS = [
    "file", "contig", "pos", "ref_in_file", "alt_in_file", "ref_in_fasta", "line_no",
]

NOT_INTERROGATED_STATEMENT = (
    "{variant} is reported in the {role} file at a position outside the interrogated "
    "regions, so the two files were not compared there. This is not a statement that the "
    "variant is present or absent in the other file."
)


def _canon(contig: str, harmonize: bool) -> str:
    if harmonize and contig.lower().startswith("chr"):
        return contig[3:]
    return contig


class _Sink:
    """Writes rows to a TSV file and/or collects them in memory."""

    def __init__(self, path: Optional[str], columns: List[str], collect: bool = False):
        self.columns = columns
        self.rows: List[Dict[str, object]] = []
        self.collect = collect
        self.n = 0
        self._fh = None
        self._writer = None
        if path:
            self._fh = open(path, "w", newline="", encoding="utf-8")
            self._writer = csv.DictWriter(self._fh, fieldnames=columns, delimiter="\t",
                                          extrasaction="ignore")
            self._writer.writeheader()

    def add(self, row: Dict[str, object]):
        self.n += 1
        if self._writer:
            self._writer.writerow(row)
        if self.collect:
            self.rows.append(row)

    def close(self):
        if self._fh:
            self._fh.close()


@dataclass
class ComparisonResult:
    """Counts, provenance and (optionally) the rows of a comparison."""

    baseline_path: str = ""
    followup_path: str = ""
    baseline_header: object = None
    followup_header: object = None
    compatibility: Optional[CompatibilityReport] = None
    baseline_stats: object = None
    followup_stats: object = None
    n_shared: int = 0
    n_baseline_only: int = 0
    n_followup_only: int = 0
    n_genotype_differences: int = 0
    baseline_only_contexts: Dict[str, int] = field(default_factory=dict)
    followup_only_contexts: Dict[str, int] = field(default_factory=dict)
    shared_rows: List[Dict[str, object]] = field(default_factory=list)
    baseline_only_rows: List[Dict[str, object]] = field(default_factory=list)
    followup_only_rows: List[Dict[str, object]] = field(default_factory=list)
    genotype_difference_rows: List[Dict[str, object]] = field(default_factory=list)
    contigs_compared: List[str] = field(default_factory=list)
    contigs_baseline_only: List[str] = field(default_factory=list)
    contigs_followup_only: List[str] = field(default_factory=list)
    normalization_level: str = NORMALIZATION_LEVEL
    options: Dict[str, object] = field(default_factory=dict)
    output_files: Dict[str, str] = field(default_factory=dict)
    # -- v0.2.0 --
    reference_path: Optional[str] = None
    reference_summary: Dict[str, object] = field(default_factory=dict)
    region_provenance: Dict[str, object] = field(default_factory=dict)
    regions_used: bool = False
    region_mode: str = "start"
    n_not_interrogated_baseline: int = 0
    n_not_interrogated_followup: int = 0
    not_interrogated_rows: List[Dict[str, object]] = field(default_factory=list)
    ref_mismatch_rows: List[Dict[str, object]] = field(default_factory=list)

    @property
    def n_reported_baseline(self) -> int:
        return self.n_shared + self.n_baseline_only

    @property
    def n_reported_followup(self) -> int:
        return self.n_shared + self.n_followup_only

    def summary(self) -> Dict[str, object]:
        extra: Dict[str, object] = {}
        if self.regions_used:
            extra["not_interrogated_baseline"] = self.n_not_interrogated_baseline
            extra["not_interrogated_followup"] = self.n_not_interrogated_followup
        if self.reference_path:
            extra["ref_mismatches"] = len(self.ref_mismatch_rows)
        return {
            "reported_in_both": self.n_shared,
            "baseline_only": self.n_baseline_only,
            "followup_only": self.n_followup_only,
            "genotype_differences_among_shared": self.n_genotype_differences,
            "reported_alleles_baseline": self.n_reported_baseline,
            "reported_alleles_followup": self.n_reported_followup,
            "jaccard": (self.n_shared / (self.n_shared + self.n_baseline_only + self.n_followup_only))
            if (self.n_shared + self.n_baseline_only + self.n_followup_only) else None,
            "normalization_level": self.normalization_level,
            **extra,
        }


def _contig_blocks(reader: VcfReader, harmonize: bool, on_not_interrogated=None
                   ) -> Iterator[Tuple[str, Dict]]:
    """Yield ``(contig, {key: AlleleCall})`` one contig at a time, in file order.

    Alleles outside the interrogated regions are handed to
    ``on_not_interrogated`` and kept out of the comparison entirely: they are
    neither shared nor unique, because the question was never asked there.
    """
    current: Optional[str] = None
    block: Dict[Tuple, AlleleCall] = {}
    finished = set()
    for call in reader:
        contig = _canon(call.contig, harmonize)
        if not call.interrogated:
            if on_not_interrogated is not None:
                on_not_interrogated(call, contig)
            continue
        if contig != current:
            if current is not None:
                if contig in finished:
                    raise UnsortedInputError(
                        f"contig {contig!r} reappears after {current!r}; records must be "
                        "grouped by contig (sort the VCF first)",
                        path=reader.path, line_no=call.line_no,
                    )
                finished.add(current)
                yield current, block
                block = {}
            current = contig
        key = (contig, call.pos, call.ref, call.alt)
        if key not in block:
            block[key] = call
    if current is not None:
        yield current, block


def _reported(call: Optional[AlleleCall], pass_only: bool, require_carried: bool) -> bool:
    """Is this allele *reported* in its file under the chosen strictness?"""
    if call is None:
        return False
    if pass_only and not call.is_pass:
        return False
    if require_carried and not call.carried:
        return False
    return True


def _context(call: Optional[AlleleCall], site_index: Dict, key, pass_only: bool,
             require_carried: bool) -> Tuple[str, str]:
    """Explain why the other file does not report this allele."""
    if call is None:
        site = (key[0], key[1])
        if site in site_index:
            return CTX_OTHER_ALLELE, ",".join(sorted(site_index[site])[:5])
        return CTX_NO_RECORD, ""
    if pass_only and not call.is_pass:
        return CTX_FILTERED, ";".join(call.filters)
    if call.gt_status == GT_NO_CALL:
        return CTX_NO_CALL, str(call.gt)
    if call.gt_status == GT_HOM_REF:
        return CTX_HOM_REF, str(call.gt)
    if not call.is_pass:
        return CTX_FILTERED, ";".join(call.filters)
    return CTX_NO_RECORD, ""


def _diff_row(call: AlleleCall, key, role: str, ctx: str, detail: str) -> Dict[str, object]:
    other = FOLLOWUP if role == BASELINE else BASELINE
    reason = CONTEXT_WORDING[ctx].format(other=other, detail=detail or "unspecified filter")
    statement = (
        f"{format_key(key)} is reported in the {role} file and is not reported in the "
        f"{other} file ({reason}). This is a difference between the two files as submitted; "
        "it is not evidence of a biological gain or loss of this variant."
    )
    return {
        "variant": format_key(key), "contig": key[0], "pos": key[1], "ref": key[2], "alt": key[3],
        "var_type": call.var_type, "reported_in": role,
        "as_written": f"{call.contig}:{call.orig_pos}:{call.orig_ref}>{call.orig_alt}",
        "id": call.var_id, "filter": ";".join(call.filters) or "PASS",
        "genotype": call.gt if call.gt is not None else "", "gt_status": call.gt_status,
        "multiallelic": call.multiallelic, "other_file_context": ctx,
        "context_detail": detail, "line_no": call.line_no, "statement": statement,
    }


def _shared_row(key, b: AlleleCall, f: AlleleCall) -> Dict[str, object]:
    """Row for an allele reported in both files.

    ``genotype_agrees`` compares the genotype *for this allele* (zygosity status
    and copy number), not the raw GT string: at a multiallelic site a GT of
    ``1/2`` and a GT of ``0/1`` can both report one copy of the same ALT.
    """
    if b.gt is None or f.gt is None:
        agrees = None
    else:
        agrees = (b.gt_status == f.gt_status and b.allele_copies == f.allele_copies)
    return {
        "variant": format_key(key), "contig": key[0], "pos": key[1], "ref": key[2], "alt": key[3],
        "var_type": b.var_type,
        "baseline_as_written": f"{b.contig}:{b.orig_pos}:{b.orig_ref}>{b.orig_alt}",
        "followup_as_written": f"{f.contig}:{f.orig_pos}:{f.orig_ref}>{f.orig_alt}",
        "baseline_filter": ";".join(b.filters) or "PASS",
        "followup_filter": ";".join(f.filters) or "PASS",
        "baseline_genotype": b.gt if b.gt is not None else "",
        "followup_genotype": f.gt if f.gt is not None else "",
        "baseline_gt_status": b.gt_status, "followup_gt_status": f.gt_status,
        "baseline_allele_copies": "" if b.allele_copies is None else b.allele_copies,
        "followup_allele_copies": "" if f.allele_copies is None else f.allele_copies,
        "genotype_agrees": "" if agrees is None else agrees,
        "baseline_line_no": b.line_no, "followup_line_no": f.line_no,
    }


def compare_files(baseline_path: str, followup_path: str, outdir: Optional[str] = None,
                  sample_baseline: Optional[str] = None, sample_followup: Optional[str] = None,
                  chr_naming: str = "strict", pass_only: bool = False,
                  require_carried: bool = False, skip_star: bool = True,
                  skip_symbolic: bool = False, emit_shared: bool = False,
                  collect: bool = False, allow_build_mismatch: bool = False,
                  require_same_sample: bool = False,
                  reference: Optional[str] = None,
                  regions: Optional[str] = None,
                  regions_baseline: Optional[str] = None,
                  regions_followup: Optional[str] = None,
                  region_mode: str = "start",
                  emit_not_interrogated: bool = False) -> ComparisonResult:
    """Compare two VCFs and return a :class:`ComparisonResult`.

    Parameters
    ----------
    outdir:
        If given, difference tables are streamed to TSV files there.
    emit_shared:
        Also write the (potentially very large) shared-variant table.
    collect:
        Keep rows in memory on the result object. Intended for tests and small
        inputs; leave off for whole-genome comparisons.
    pass_only:
        Count only FILTER=PASS records as reported.
    require_carried:
        Count an allele as reported only when the sample genotype carries it.
    """
    harmonize = chr_naming == "harmonize"
    fasta = IndexedFasta(reference) if reference else None
    region_set, region_prov = build_comparison_regions(regions, regions_baseline,
                                                       regions_followup)
    b_reader = VcfReader(baseline_path, sample=sample_baseline, skip_star=skip_star,
                         skip_symbolic=skip_symbolic, fasta=fasta, regions=region_set,
                         region_mode=region_mode)
    f_reader = VcfReader(followup_path, sample=sample_followup, skip_star=skip_star,
                         skip_symbolic=skip_symbolic, fasta=fasta, regions=region_set,
                         region_mode=region_mode)
    compat = check_compatibility(b_reader.header, f_reader.header, chr_naming=chr_naming,
                                require_same_sample=require_same_sample,
                                allow_build_mismatch=allow_build_mismatch)
    try:
        require_compatible(compat)
    except Exception:
        b_reader.close()
        f_reader.close()
        raise

    res = ComparisonResult(
        baseline_path=baseline_path, followup_path=followup_path,
        baseline_header=b_reader.header, followup_header=f_reader.header,
        compatibility=compat,
        options={"chr_naming": chr_naming, "pass_only": pass_only,
                 "require_carried": require_carried, "skip_star": skip_star,
                 "skip_symbolic": skip_symbolic, "emit_shared": emit_shared,
                 "sample_baseline": b_reader.sample, "sample_followup": f_reader.sample},
        normalization_level=LEVEL_LEFT_ALIGNED if fasta else LEVEL_TRIM_ONLY,
        reference_path=reference,
        reference_summary=fasta.summary() if fasta else {},
        region_provenance=region_prov,
        regions_used=region_set is not None,
        region_mode=region_mode,
    )
    if fasta:
        res.options["reference"] = os.path.basename(reference)
    if region_set is not None:
        res.options["regions"] = region_prov.get("combination", "")
        res.options["region_mode"] = region_mode
    sinks = {}
    if outdir:
        os.makedirs(outdir, exist_ok=True)
    def _p(name):
        return os.path.join(outdir, name) if outdir else None
    sinks["baseline_only"] = _Sink(_p("baseline_only.tsv"), DIFF_COLUMNS, collect)
    sinks["followup_only"] = _Sink(_p("followup_only.tsv"), DIFF_COLUMNS, collect)
    sinks["shared"] = _Sink(_p("shared.tsv") if emit_shared else None, SHARED_COLUMNS, collect)
    sinks["genotype_differences"] = _Sink(_p("genotype_differences.tsv"), SHARED_COLUMNS, collect)
    want_ni = region_set is not None and emit_not_interrogated
    sinks["not_interrogated"] = _Sink(_p("not_interrogated.tsv") if want_ni else None,
                                      NOT_INTERROGATED_COLUMNS, collect and region_set is not None)

    def _note_not_interrogated(call: AlleleCall, contig: str, role: str):
        if role == BASELINE:
            res.n_not_interrogated_baseline += 1
        else:
            res.n_not_interrogated_followup += 1
        if sinks["not_interrogated"].collect or want_ni:
            key = (contig, call.pos, call.ref, call.alt)
            sinks["not_interrogated"].add({
                "variant": format_key(key), "contig": contig, "pos": call.pos,
                "ref": call.ref, "alt": call.alt, "var_type": call.var_type,
                "reported_in": role,
                "as_written": f"{call.contig}:{call.orig_pos}:{call.orig_ref}>{call.orig_alt}",
                "filter": ";".join(call.filters) or "PASS",
                "genotype": call.gt if call.gt is not None else "",
                "line_no": call.line_no,
                "statement": NOT_INTERROGATED_STATEMENT.format(
                    variant=format_key(key), role=role),
            })

    def _emit_pair(b_block: Dict, f_block: Dict):
        b_sites: Dict[Tuple, List[str]] = {}
        f_sites: Dict[Tuple, List[str]] = {}
        for k in b_block:
            b_sites.setdefault((k[0], k[1]), []).append(f"{k[2]}>{k[3]}")
        for k in f_block:
            f_sites.setdefault((k[0], k[1]), []).append(f"{k[2]}>{k[3]}")
        for key, b in b_block.items():
            f = f_block.get(key)
            b_rep = _reported(b, pass_only, require_carried)
            f_rep = _reported(f, pass_only, require_carried)
            if b_rep and f_rep:
                row = _shared_row(key, b, f)
                res.n_shared += 1
                sinks["shared"].add(row)
                if row["genotype_agrees"] is False:
                    res.n_genotype_differences += 1
                    sinks["genotype_differences"].add(row)
            elif b_rep and not f_rep:
                ctx, detail = _context(f, f_sites, key, pass_only, require_carried)
                res.n_baseline_only += 1
                res.baseline_only_contexts[ctx] = res.baseline_only_contexts.get(ctx, 0) + 1
                sinks["baseline_only"].add(_diff_row(b, key, BASELINE, ctx, detail))
            elif f_rep and not b_rep:
                ctx, detail = _context(b, b_sites, key, pass_only, require_carried)
                res.n_followup_only += 1
                res.followup_only_contexts[ctx] = res.followup_only_contexts.get(ctx, 0) + 1
                sinks["followup_only"].add(_diff_row(f, key, FOLLOWUP, ctx, detail))
        for key, f in f_block.items():
            if key in b_block:
                continue
            if not _reported(f, pass_only, require_carried):
                continue
            ctx, detail = _context(None, b_sites, key, pass_only, require_carried)
            res.n_followup_only += 1
            res.followup_only_contexts[ctx] = res.followup_only_contexts.get(ctx, 0) + 1
            sinks["followup_only"].add(_diff_row(f, key, FOLLOWUP, ctx, detail))

    try:
        b_blocks = _contig_blocks(
            b_reader, harmonize,
            on_not_interrogated=(lambda c, k: _note_not_interrogated(c, k, BASELINE))
            if region_set is not None else None)
        f_blocks = _contig_blocks(
            f_reader, harmonize,
            on_not_interrogated=(lambda c, k: _note_not_interrogated(c, k, FOLLOWUP))
            if region_set is not None else None)
        b_pending: Dict[str, Dict] = {}
        f_pending: Dict[str, Dict] = {}
        b_done = f_done = False
        while not (b_done and f_done):
            if not b_done:
                try:
                    contig, block = next(b_blocks)
                    if contig in f_pending:
                        res.contigs_compared.append(contig)
                        _emit_pair(block, f_pending.pop(contig))
                    else:
                        b_pending[contig] = block
                except StopIteration:
                    b_done = True
            if not f_done:
                try:
                    contig, block = next(f_blocks)
                    if contig in b_pending:
                        res.contigs_compared.append(contig)
                        _emit_pair(b_pending.pop(contig), block)
                    else:
                        f_pending[contig] = block
                except StopIteration:
                    f_done = True
        for contig, block in b_pending.items():
            res.contigs_baseline_only.append(contig)
            _emit_pair(block, {})
        for contig, block in f_pending.items():
            res.contigs_followup_only.append(contig)
            _emit_pair({}, block)
    finally:
        b_reader.close()
        f_reader.close()
        for name, sink in sinks.items():
            sink.close()
            skip = (name == "shared" and not emit_shared) or \
                   (name == "not_interrogated" and not want_ni)
            if outdir and not skip:
                res.output_files[name] = _p(f"{name}.tsv")

    res.baseline_stats = b_reader.stats
    res.followup_stats = f_reader.stats

    # REF-vs-reference disagreements are written out, not hidden in a counter.
    if fasta is not None:
        for role, stats in ((BASELINE, b_reader.stats), (FOLLOWUP, f_reader.stats)):
            for ex in stats.ref_mismatch_examples:
                res.ref_mismatch_rows.append(dict(file=role, **ex))
        if res.ref_mismatch_rows and outdir:
            path = _p("ref_mismatches.tsv")
            with open(path, "w", newline="", encoding="utf-8") as fh:
                w = csv.DictWriter(fh, fieldnames=REF_MISMATCH_COLUMNS, delimiter="\t",
                                   extrasaction="ignore")
                w.writeheader()
                for row in res.ref_mismatch_rows:
                    w.writerow(row)
            res.output_files["ref_mismatches"] = path
        fasta.close()

    if collect:
        res.shared_rows = sinks["shared"].rows
        res.baseline_only_rows = sinks["baseline_only"].rows
        res.followup_only_rows = sinks["followup_only"].rows
        res.genotype_difference_rows = sinks["genotype_differences"].rows
        res.not_interrogated_rows = sinks["not_interrogated"].rows
    return res
