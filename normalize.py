"""Allele normalization.

GTT has two normalization levels.

``trim_only`` (no reference needed)
    *Parsimonious trimming*: common suffix then common prefix are removed,
    leaving at least one base in each allele. Deterministic, but NOT full
    left-alignment: left-shifting an indel through a repeat requires the
    reference sequence. Two files that anchor the same indel at different
    positions inside a homopolymer or tandem repeat appear as one baseline-only
    and one follow-up-only variant.

``left_aligned`` (``--reference`` given)
    Full left-alignment by the algorithm of Tan, Abecasis & Kang (2015),
    *Unified representation of genetic variants*, Bioinformatics 31(13):2202-4 -
    the algorithm behind ``vt normalize`` and ``bcftools norm``. An indel is
    slid as far left as the reference sequence allows, so that every equivalent
    way of writing it collapses onto one representation.

The algorithm in plain language
-------------------------------
Take the deletion of one A from the run ``C AAAAAAAA G``. It can be written as
"drop the A at position 11", "drop the A at position 17", and six other ways:
all describe the same resulting sequence. Left-alignment picks the leftmost one
by repeating two moves:

1. if both alleles end with the same base, that base carries no information, so
   remove it from both (first lengthening the alleles one base to the left when
   one of them is down to a single base, so nothing becomes empty);
2. stop when the last bases differ, then drop any shared leading bases.

Running this on every file makes equivalent indels identical, which is the
whole point: a difference table should list real disagreements, not two
spellings of the same thing.
"""
from __future__ import annotations

import re
from typing import Optional, Tuple

NORMALIZATION_LEVEL = "trim_only"
LEVEL_TRIM_ONLY = "trim_only"
LEVEL_LEFT_ALIGNED = "left_aligned"

#: Furthest an allele is shifted left before GTT gives up and says so.
MAX_LEFT_SHIFT = 10000

# Outcome of a left-alignment attempt, recorded per allele and counted.
LA_UNCHANGED = "unchanged"                      # already left-aligned
LA_SHIFTED = "shifted"                          # moved and/or re-spelled
LA_NOT_APPLICABLE = "not_applicable"            # symbolic, '*', '.' - nothing to align
LA_REF_MISMATCH = "ref_mismatch"                # file REF != reference FASTA; left as trimmed
LA_CONTIG_ABSENT = "contig_not_in_reference"    # contig missing from the FASTA
LA_AT_CONTIG_START = "truncated_at_contig_start"  # cannot shift past base 1
LA_SHIFT_LIMIT = "shift_limit_reached"          # hit MAX_LEFT_SHIFT

#: Statuses that mean "the reference could not be used for this allele".
LA_PROBLEM_STATUSES = (LA_REF_MISMATCH, LA_CONTIG_ABSENT, LA_SHIFT_LIMIT)

_DNA_RE = re.compile(r"^[ACGTN]+$")
_SYMBOLIC_RE = re.compile(r"^<[^<>]+>$")
_BREAKEND_RE = re.compile(r"[\[\]]")

#: ALT values that do not denote an independent sequence change.
NON_VARIANT_ALTS = {".", "*"}


def is_symbolic(alt: str) -> bool:
    """True for ``<DEL>``-style symbolic alleles and breakend notation."""
    return bool(_SYMBOLIC_RE.match(alt)) or bool(_BREAKEND_RE.search(alt))


def is_plain_dna(allele: str) -> bool:
    """True when the allele is a plain [ACGTN] sequence."""
    return bool(_DNA_RE.match(allele))


def normalize_allele(pos: int, ref: str, alt: str) -> Tuple[int, str, str]:
    """Return the parsimoniously trimmed ``(pos, ref, alt)``.

    Symbolic alleles, breakends and the ``*`` / ``.`` placeholders are returned
    uppercased but otherwise untouched, since trimming them is not meaningful.

    >>> normalize_allele(100, "CTT", "CT")       # 1 bp deletion
    (101, 'TT', 'T')
    >>> normalize_allele(100, "T", "TGCCCA")     # insertion, already minimal
    (100, 'T', 'TGCCCA')
    >>> normalize_allele(100, "AGT", "CGT")      # trims to an SNV
    (100, 'A', 'C')
    """
    ref = ref.upper()
    alt = alt.upper()
    if alt in NON_VARIANT_ALTS or is_symbolic(alt) or not is_plain_dna(ref):
        return pos, ref, alt
    if not is_plain_dna(alt):
        return pos, ref, alt

    # Trim the common suffix first (leave >= 1 base in each allele).
    while len(ref) > 1 and len(alt) > 1 and ref[-1] == alt[-1]:
        ref = ref[:-1]
        alt = alt[:-1]
    # Then trim the common prefix, advancing POS for every base removed.
    trim = 0
    while trim < len(ref) - 1 and trim < len(alt) - 1 and ref[trim] == alt[trim]:
        trim += 1
    if trim:
        ref = ref[trim:]
        alt = alt[trim:]
        pos += trim
    return pos, ref, alt


def left_align_allele(fasta, contig: str, pos: int, ref: str, alt: str,
                      max_shift: int = MAX_LEFT_SHIFT, window: int = 64
                      ) -> Tuple[int, str, str, str, Optional[str]]:
    """Left-align one allele against a reference FASTA.

    ``pos``/``ref``/``alt`` should already be parsimoniously trimmed (the
    output of :func:`normalize_allele`). ``fasta`` is an
    :class:`~gtt.fasta.IndexedFasta`.

    Returns ``(pos, ref, alt, status, ref_in_fasta)``. ``status`` is one of the
    ``LA_*`` constants. ``ref_in_fasta`` is filled only for ``ref_mismatch``
    and holds what the reference actually says at that position, so the
    mismatch can be reported rather than hidden.

    The REF allele is checked against the reference first. A mismatch means the
    VCF and the FASTA are not describing the same genome at that position
    (wrong build, wrong patch, or a corrupted file); the allele is then left
    exactly as trimmed and counted, never silently "corrected".
    """
    if alt in NON_VARIANT_ALTS or is_symbolic(alt) or not is_plain_dna(ref) \
            or not is_plain_dna(alt):
        return pos, ref, alt, LA_NOT_APPLICABLE, None
    if ref == alt:
        return pos, ref, alt, LA_NOT_APPLICABLE, None

    name = fasta.resolve(contig)
    if name is None:
        return pos, ref, alt, LA_CONTIG_ABSENT, None
    contig_length = fasta.records[name].length

    # 1. Does the file's REF match the reference sequence?
    if pos < 1 or pos - 1 + len(ref) > contig_length:
        return pos, ref, alt, LA_REF_MISMATCH, None
    observed = fasta.fetch(name, pos - 1, pos - 1 + len(ref))
    if observed != ref:
        return pos, ref, alt, LA_REF_MISMATCH, observed

    # 2. Tan/Abecasis/Kang left-shift.
    start_pos = pos
    original = (pos, ref, alt)
    status = LA_UNCHANGED
    buf = ""          # reference bases immediately left of `pos`
    while True:
        if ref[-1] != alt[-1]:
            break
        if len(ref) == 1 or len(alt) == 1:
            if pos == 1:
                status = LA_AT_CONTIG_START
                break
            if start_pos - pos >= max_shift:
                status = LA_SHIFT_LIMIT
                break
            if not buf:
                w = min(window, pos - 1)
                buf = fasta.fetch(name, pos - 1 - w, pos - 1)
            base, buf = buf[-1], buf[:-1]
            pos -= 1
            ref = base + ref
            alt = base + alt
        ref = ref[:-1]
        alt = alt[:-1]
    # 3. Drop shared leading bases (keep at least one base in each allele).
    while len(ref) >= 2 and len(alt) >= 2 and ref[0] == alt[0]:
        ref = ref[1:]
        alt = alt[1:]
        pos += 1
    if status == LA_UNCHANGED and (pos, ref, alt) != original:
        status = LA_SHIFTED
    return pos, ref, alt, status, None


def classify_variant(ref: str, alt: str) -> str:
    """Classify a normalized allele pair.

    Returns one of ``snv``, ``mnv``, ``insertion``, ``deletion``, ``indel``,
    ``symbolic``, ``upstream_deletion`` (``*``) or ``no_alt`` (``.``).
    """
    if alt == "*":
        return "upstream_deletion"
    if alt == ".":
        return "no_alt"
    if is_symbolic(alt):
        return "symbolic"
    if len(ref) == len(alt):
        return "snv" if len(ref) == 1 else "mnv"
    if len(ref) == 1 and alt.startswith(ref):
        return "insertion"
    if len(alt) == 1 and ref.startswith(alt):
        return "deletion"
    return "indel"
