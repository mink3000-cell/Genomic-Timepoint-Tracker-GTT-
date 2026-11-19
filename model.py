"""Core data model: one row per (record, ALT allele) after decomposition."""
from __future__ import annotations

from dataclasses import dataclass
from typing import Optional, Tuple

# Genotype status vocabulary (per ALT allele, for one sample).
GT_NO_SAMPLE = "no_sample"      # VCF has no sample columns at all
GT_NO_CALL = "no_call"          # GT present but missing (./.)
GT_HOM_REF = "hom_ref"          # GT is 0/0 - reference call, allele not carried
GT_HET = "het"                  # one copy of this ALT
GT_HOM_ALT = "hom_alt"          # all called alleles are this ALT
GT_HEMI = "hemizygous"          # haploid GT carrying this ALT
GT_OTHER_ALT = "other_alt_only" # GT carries a different ALT at a multiallelic site
GT_NO_GT_FIELD = "no_gt_field"  # FORMAT has no GT key


@dataclass(slots=True)
class AlleleCall:
    """A single ALT allele of a single VCF record, after decomposition.

    Coordinates in ``contig``/``pos``/``ref``/``alt`` are normalized; the
    ``orig_*`` fields preserve the representation as written in the file so the
    report can always be traced back to its source line.
    """

    contig: str
    pos: int
    ref: str
    alt: str
    orig_pos: int
    orig_ref: str
    orig_alt: str
    var_id: str = "."
    qual: str = "."
    filters: Tuple[str, ...] = ()
    is_pass: bool = True
    gt: Optional[str] = None
    gt_status: str = GT_NO_SAMPLE
    allele_copies: Optional[int] = None
    phased: bool = False
    multiallelic: bool = False
    n_alts: int = 1
    alt_index: int = 1
    var_type: str = "snv"
    symbolic: bool = False
    line_no: int = 0
    sample: Optional[str] = None
    #: False when the allele falls outside the interrogated regions (``--regions``).
    #: Such alleles are reported as *not interrogated*, never as absent.
    interrogated: bool = True
    #: Outcome of left-alignment for this allele (``""`` when no reference was used).
    la_status: str = ""

    @property
    def key(self):
        """Comparison key: normalized contig/pos/ref/alt. Genotype is deliberately excluded."""
        return (self.contig, self.pos, self.ref, self.alt)

    @property
    def site(self):
        """Site key: normalized contig/pos only."""
        return (self.contig, self.pos)

    @property
    def carried(self) -> bool:
        """True when the sample genotype carries at least one copy of this ALT."""
        return self.gt_status in (GT_HET, GT_HOM_ALT, GT_HEMI)


def variant_key(contig, pos, ref, alt):
    """Build a comparison key from already-normalized parts."""
    return (str(contig), int(pos), str(ref).upper(), str(alt).upper())


def format_key(key) -> str:
    """Human-readable rendering of a comparison key, e.g. ``chr1-604358-G-C``."""
    contig, pos, ref, alt = key
    return f"{contig}-{pos}-{ref}-{alt}"
