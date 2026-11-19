"""Input validation and reference-genome compatibility checks.

Cross-build comparison is refused rather than warned about: comparing GRCh37
coordinates against GRCh38 coordinates yields a report that is almost entirely
spurious differences. GTT has no liftover, so it will not pretend to.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from .errors import BuildMismatchError, VcfFormatError
from .vcfio import VcfHeader

#: Length of chromosome 1 in each supported assembly - a unique fingerprint.
_CHR1_LENGTHS = {
    248956422: "GRCh38",
    249250621: "GRCh37",
    247249719: "NCBI36",
    248387328: "T2T-CHM13v2.0",
}
#: Secondary fingerprint on chromosome 2, used when chromosome 1 is absent.
_CHR2_LENGTHS = {
    242193529: "GRCh38",
    243199373: "GRCh37",
    242951149: "NCBI36",
    242696752: "T2T-CHM13v2.0",
}

_ASSEMBLY_HINTS = (
    ("grch38", "GRCh38"), ("hg38", "GRCh38"),
    ("grch37", "GRCh37"), ("hg19", "GRCh37"), ("b37", "GRCh37"),
    ("chm13", "T2T-CHM13v2.0"), ("hg18", "NCBI36"), ("ncbi36", "NCBI36"),
)


def _canonical(name: str) -> str:
    """Strip a ``chr`` prefix so ``chr1`` and ``1`` compare equal."""
    return name[3:] if name.lower().startswith("chr") else name


def infer_build(header: VcfHeader):
    """Infer the assembly from contig-length fingerprints, then from name hints.

    Returns ``(build, evidence)``; ``build`` is None when nothing is decisive.
    """
    by_canon = {_canonical(c): length for c, length in header.contigs.items()}
    for chrom, table in (("1", _CHR1_LENGTHS), ("2", _CHR2_LENGTHS)):
        length = by_canon.get(chrom)
        if length in table:
            return table[length], f"contig {chrom} length {length}"
    hint_sources = [header.reference or ""] + list(header.contig_assembly.values())
    for src in hint_sources:
        low = src.lower()
        for needle, build in _ASSEMBLY_HINTS:
            if needle in low:
                return build, f"assembly/reference string {src!r}"
    return None, "no contig-length fingerprint and no recognizable assembly string"


@dataclass
class CompatibilityReport:
    """Result of comparing two headers for comparability."""

    baseline_build: Optional[str] = None
    followup_build: Optional[str] = None
    baseline_evidence: str = ""
    followup_evidence: str = ""
    chr_naming_mismatch: bool = False
    length_conflicts: List[str] = field(default_factory=list)
    shared_contigs: List[str] = field(default_factory=list)
    baseline_only_contigs: List[str] = field(default_factory=list)
    followup_only_contigs: List[str] = field(default_factory=list)
    sample_note: str = ""
    warnings: List[str] = field(default_factory=list)
    compatible: bool = True
    refusal_reason: Optional[str] = None

    def as_dict(self) -> Dict[str, object]:
        return {k: getattr(self, k) for k in self.__dataclass_fields__}


def validate_header(header: VcfHeader) -> List[str]:
    """Structural checks on a parsed header. Returns a list of warnings."""
    warnings: List[str] = []
    if header.fileformat and not header.fileformat.startswith("VCFv4"):
        warnings.append(
            f"{header.path}: ##fileformat is {header.fileformat!r}; GTT is tested against VCFv4.x"
        )
    if not header.contigs:
        warnings.append(
            f"{header.path}: no ##contig lines - reference build cannot be fingerprinted "
            "from contig lengths"
        )
    if header.n_samples == 0:
        warnings.append(f"{header.path}: no sample columns - genotypes will not be compared")
    elif header.n_samples > 1:
        warnings.append(
            f"{header.path}: {header.n_samples} samples present ({', '.join(header.samples)}); "
            "the first is used unless --sample is given"
        )
    if header.synthetic:
        warnings.append(
            f"{header.path}: SYNTHETIC TEST FIXTURE ({header.synthetic_note}) - "
            "not a biological observation"
        )
    return warnings


def check_compatibility(baseline: VcfHeader, followup: VcfHeader,
                        chr_naming: str = "strict", require_same_sample: bool = False,
                        allow_build_mismatch: bool = False) -> CompatibilityReport:
    """Decide whether two VCFs may be compared.

    ``chr_naming`` is ``strict`` (refuse a ``chr1`` vs ``1`` mismatch) or
    ``harmonize`` (compare on canonical names, with a warning).
    """
    rep = CompatibilityReport()
    rep.baseline_build, rep.baseline_evidence = infer_build(baseline)
    rep.followup_build, rep.followup_evidence = infer_build(followup)
    rep.warnings.extend(validate_header(baseline))
    rep.warnings.extend(validate_header(followup))

    if rep.baseline_build and rep.followup_build and rep.baseline_build != rep.followup_build:
        rep.compatible = False
        rep.refusal_reason = (
            f"reference build mismatch: baseline is {rep.baseline_build} "
            f"({rep.baseline_evidence}), follow-up is {rep.followup_build} "
            f"({rep.followup_evidence}). GTT does not lift over coordinates."
        )
    elif not rep.baseline_build or not rep.followup_build:
        rep.warnings.append(
            "reference build could not be confirmed for at least one input "
            f"(baseline: {rep.baseline_evidence}; follow-up: {rep.followup_evidence}); "
            "coordinate comparability is assumed, not verified"
        )

    b_pref, f_pref = baseline.chr_prefixed, followup.chr_prefixed
    if b_pref is not None and f_pref is not None and b_pref != f_pref:
        rep.chr_naming_mismatch = True
        if chr_naming == "strict":
            rep.compatible = False
            rep.refusal_reason = (
                (rep.refusal_reason + " | " if rep.refusal_reason else "")
                + "contig naming mismatch ('chr1' vs '1'). Re-run with "
                  "--chr-naming harmonize to compare on canonical names."
            )
        else:
            rep.warnings.append(
                "contig naming differs between inputs; comparing on canonical names "
                "(chr prefix ignored)"
            )

    b_map = {_canonical(c): l for c, l in baseline.contigs.items()}
    f_map = {_canonical(c): l for c, l in followup.contigs.items()}
    rep.shared_contigs = [c for c in b_map if c in f_map]
    rep.baseline_only_contigs = [c for c in b_map if c not in f_map]
    rep.followup_only_contigs = [c for c in f_map if c not in b_map]
    for c in rep.shared_contigs:
        if b_map[c] is not None and f_map[c] is not None and b_map[c] != f_map[c]:
            rep.length_conflicts.append(f"{c}: baseline {b_map[c]} vs follow-up {f_map[c]}")
    if rep.length_conflicts:
        rep.compatible = False
        rep.refusal_reason = (
            (rep.refusal_reason + " | " if rep.refusal_reason else "")
            + "contig length conflicts: " + "; ".join(rep.length_conflicts[:5])
        )

    if baseline.samples and followup.samples:
        if baseline.samples[0] != followup.samples[0]:
            rep.sample_note = (
                f"sample identifiers differ (baseline {baseline.samples[0]!r}, "
                f"follow-up {followup.samples[0]!r}); GTT cannot verify that the two files "
                "describe the same individual"
            )
            if require_same_sample:
                rep.compatible = False
                rep.refusal_reason = (
                    (rep.refusal_reason + " | " if rep.refusal_reason else "")
                    + rep.sample_note
                )
            else:
                rep.warnings.append(rep.sample_note)
        else:
            rep.sample_note = f"both files carry sample {baseline.samples[0]!r}"

    if not rep.compatible and allow_build_mismatch:
        rep.warnings.append(
            "OVERRIDE: --allow-build-mismatch given; proceeding despite: " + str(rep.refusal_reason)
        )
        rep.compatible = True
    return rep


def require_compatible(rep: CompatibilityReport) -> None:
    """Raise :class:`BuildMismatchError` if the inputs are not comparable."""
    if not rep.compatible:
        raise BuildMismatchError(rep.refusal_reason or "inputs are not comparable")
