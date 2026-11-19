"""Report rendering with an enforced, conservative vocabulary.

Wording rules, applied in code rather than left to the writer:

* a variant missing from a file is "not reported in <file>" - never "lost",
  "reverted", "disappeared" or "cleared";
* a variant seen only in the follow-up file is "reported only in the follow-up
  file" - never "newly acquired", "de novo" or "treatment-induced";
* every difference count is accompanied by the technical reasons that can
  produce it (coverage, filtering, caller version, normalization).
"""
from __future__ import annotations

import json
import os
from datetime import datetime, timezone
from typing import Dict, Optional

from . import __version__
from .compare import (
    CTX_FILTERED,
    CTX_HOM_REF,
    CTX_NO_CALL,
    CTX_NO_RECORD,
    CTX_OTHER_ALLELE,
    ComparisonResult,
)

#: Words that must never appear in a GTT report when describing absence.
FORBIDDEN_INTERPRETIVE_TERMS = (
    "lost", "loss of variant", "disappeared", "reverted", "cured", "cleared",
    "newly acquired", "de novo", "treatment-induced", "gained a mutation",
)

CAVEAT_TEXT = """A difference between two VCF files is a difference between two
*reports*, not an observed biological change. Any of the following can place a
variant in one file and not the other without any change in the underlying
genome: unequal sequencing depth or coverage at the position; different variant
callers, caller versions or filtering thresholds; different confident-region
definitions, so that a position is interrogated in one run and not the other;
and different but equivalent representations of the same indel. GTT reports the
file contents and the reason the other file does not report each variant; it
does not and cannot establish that a variant was gained or lost."""

NORMALIZATION_CAVEAT = """Alleles were normalized by parsimonious trimming only
(common suffix, then common prefix, POS adjusted). Full left-alignment requires
a reference FASTA, which was not supplied. Consequence: the same indel anchored
at different positions inside a homopolymer or tandem repeat will appear once in
each of the two difference tables rather than as a shared variant. Indel-heavy
difference counts should be read with that in mind."""

LEFT_ALIGNED_CAVEAT = """Alleles were fully left-aligned against the supplied
reference FASTA (parsimonious trimming, then the left-shifting algorithm of Tan,
Abecasis & Kang 2015, the algorithm used by vt normalize and bcftools norm).
Equivalent spellings of the same indel inside a homopolymer or tandem repeat
therefore collapse onto one representation and are counted as shared rather than
as one difference in each table. Two limits remain: an indel whose leftmost
representation would start before base 1 of its contig is reported at base 1 and
counted as truncated at the contig start, and a REF allele that does not match
the reference sequence is left exactly as the file wrote it, counted, and listed
in the mismatch table."""

REGION_CAVEAT = """The comparison was restricted to the supplied regions. A
variant reported outside them is listed as *not interrogated*: the two files
were not compared at that position, which is a statement about the scope of this
run and not about the presence or absence of the variant in either file. When
two region files are supplied, the comparison region is their intersection, so
that every compared position is one both files make a statement about."""

CONTEXT_LABELS = {
    CTX_NO_RECORD: "no record at that position in the other file",
    CTX_OTHER_ALLELE: "a record exists at that position but reports a different allele",
    CTX_FILTERED: "the allele is present in the other file but its record is filtered",
    CTX_NO_CALL: "the allele is present in the other file but the genotype is a no-call",
    CTX_HOM_REF: "the allele is present in the other file but the genotype is homozygous reference",
}


def _fmt(n) -> str:
    return f"{n:,}" if isinstance(n, int) else str(n)


def _header_block(role: str, header, stats) -> str:
    s = header.summary()
    lines = [
        f"**{role}** - `{s['path']}`",
        "",
        f"| field | value |", "| --- | --- |",
        f"| ##fileformat | {s['fileformat']} |",
        f"| ##reference | {s['reference'] or 'not declared'} |",
        f"| contigs declared | {_fmt(s['n_contigs'])} |",
        f"| sample column used | {', '.join(s['samples']) or 'none'} |",
        f"| records read | {_fmt(stats.records)} |",
        f"| ALT alleles after decomposition | {_fmt(stats.alleles_emitted)} |",
        f"| multiallelic records | {_fmt(stats.multiallelic_records)} |",
        f"| records with FILTER other than PASS | {_fmt(stats.not_pass_records)} |",
        f"| records with a no-call genotype | {_fmt(stats.no_call_records)} |",
        f"| `*` alleles skipped | {_fmt(stats.skipped_star)} |",
    ]
    if header.synthetic:
        note = (header.synthetic_note or "").rstrip(". ")
        suffix = "" if "not a biological observation" in note.lower() \
            else ". Not a biological observation"
        lines.insert(1, "")
        lines.insert(2, f"> SYNTHETIC TEST FIXTURE - {note}{suffix}.")
    return "\n".join(lines)


def _context_table(counts: Dict[str, int], other: str) -> str:
    if not counts:
        return "_none_"
    rows = ["| reason the variant is not reported in the " + other + " file | variants |",
            "| --- | ---: |"]
    for ctx, n in sorted(counts.items(), key=lambda kv: -kv[1]):
        rows.append(f"| {CONTEXT_LABELS.get(ctx, ctx)} | {_fmt(n)} |")
    return "\n".join(rows)


def _regions_section(res: ComparisonResult):
    """Section 4.3 - what was interrogated, and what was therefore not compared."""
    prov = res.region_provenance
    lines = ["### 4.3 Region restriction (what was interrogated)", "",
             "| region file | contigs | intervals | bases |", "| --- | ---: | ---: | ---: |"]
    for key in ("regions", "regions_baseline", "regions_followup", "comparison_regions"):
        s = prov.get(key)
        if not s:
            continue
        label = "**comparison region**" if key == "comparison_regions" else key
        src = ", ".join(s["sources"]) if s.get("sources") else ""
        lines.append(f"| {label} `{src}` | {_fmt(s['n_contigs'])} | "
                     f"{_fmt(s['n_intervals'])} | {_fmt(s['total_bases'])} |")
    nb, nf = res.n_not_interrogated_baseline, res.n_not_interrogated_followup
    lines += [
        "",
        f"Position test: `{res.region_mode}` "
        + ("(the normalized POS must be inside an interval)" if res.region_mode == "start"
           else "(every base of REF must be inside one interval)"),
        "",
        "| file | alleles not interrogated |", "| --- | ---: |",
        f"| baseline | {_fmt(nb)} |",
        f"| follow-up | {_fmt(nf)} |",
        "",
        "These alleles fall outside the interrogated region. They were not compared and are "
        "not counted in section 4; listing them as not interrogated is a statement about this "
        "comparison, not about the variant.",
        "",
    ]
    for role, stats in (("baseline", res.baseline_stats), ("follow-up", res.followup_stats)):
        by_type = getattr(stats, "not_interrogated_by_type", {}) or {}
        if by_type:
            parts = ", ".join(f"{t}: {_fmt(n)}" for t, n in sorted(by_type.items()))
            lines.append(f"- not interrogated in the {role} file by type - {parts}")
    lines.append("")
    return lines


def _reference_section(res: ComparisonResult):
    """Section 4.4 - what the reference FASTA changed, and what it flagged."""
    rs = res.reference_summary or {}
    lines = ["### 4.4 Reference check and left-alignment", "",
             "| field | value |", "| --- | --- |",
             f"| reference FASTA | `{rs.get('path', res.reference_path)}` |",
             f"| index | `{rs.get('fai', '')}`"
             f"{' (created by GTT)' if rs.get('index_created_by_gtt') else ''} |",
             f"| contigs in reference | {_fmt(rs.get('n_contigs', 0))} |",
             f"| normalization level | `{res.normalization_level}` |",
             "",
             "| outcome per ALT allele | baseline | follow-up |",
             "| --- | ---: | ---: |"]
    rows = [
        ("already left-aligned", "left_aligned_unchanged"),
        ("moved by left-alignment", "left_aligned_shifted"),
        ("REF disagrees with the reference FASTA", "ref_mismatches"),
        ("reported at base 1, could not shift further left", "truncated_at_contig_start"),
        ("shift limit reached", "shift_limit_reached"),
        ("contig absent from the reference FASTA", "alleles_on_contigs_not_in_reference"),
    ]
    for label, attr in rows:
        b = getattr(res.baseline_stats, attr, 0)
        f = getattr(res.followup_stats, attr, 0)
        lines.append(f"| {label} | {_fmt(b)} | {_fmt(f)} |")
    lines.append("")
    missing = sorted(set(getattr(res.baseline_stats, "contigs_not_in_reference", []))
                     | set(getattr(res.followup_stats, "contigs_not_in_reference", [])))
    if missing:
        lines += [
            "**Contigs absent from the reference FASTA**: " + ", ".join(missing)
            + ". Alleles on them were trimmed but not left-aligned, so the "
              "`trim_only` limitation still applies to those positions.", "",
        ]
    if res.ref_mismatch_rows:
        lines += [
            f"**{_fmt(len(res.ref_mismatch_rows))} REF "
            f"{'allele disagrees' if len(res.ref_mismatch_rows) == 1 else 'alleles disagree'} "
            "with the reference sequence.** A disagreement means the VCF and the FASTA are "
            "not describing the "
            "same genome at that position - a different assembly or patch, or a damaged file. "
            "GTT leaves those alleles exactly as the file wrote them, counts them here, and "
            "does not left-align them. First examples:", "",
            "| file | position | REF in file | REF in FASTA |", "| --- | --- | --- | --- |",
        ]
        for row in res.ref_mismatch_rows[:10]:
            lines.append(f"| {row['file']} | {row['contig']}:{row['pos']} | "
                         f"`{row['ref_in_file']}` | `{row['ref_in_fasta']}` |")
        lines.append("")
    return lines


def render_markdown(res: ComparisonResult, title: str = "GTT comparison report",
                    artifact_links: Optional[Dict[str, str]] = None) -> str:
    """Render the full comparison report as Markdown."""
    c = res.compatibility
    synthetic = bool(getattr(res.baseline_header, "synthetic", False)
                     or getattr(res.followup_header, "synthetic", False))
    out = [f"# {title}", ""]
    if synthetic:
        out += ["> **Synthetic test data.** At least one input is a GTT test fixture carrying a",
                "> `##gtt_synthetic` header. The comparison below exercises the comparison logic.",
                "> It is not a biological observation and must not be reported as one.", ""]
    out += [
        f"_Generated {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')} by GTT v{__version__}._",
        "",
        "## 1. Inputs", "",
        _header_block("Baseline", res.baseline_header, res.baseline_stats), "",
        _header_block("Follow-up", res.followup_header, res.followup_stats), "",
        "## 2. Reference genome compatibility", "",
        f"| check | result |", "| --- | --- |",
        f"| baseline build | {c.baseline_build or 'undetermined'} ({c.baseline_evidence}) |",
        f"| follow-up build | {c.followup_build or 'undetermined'} ({c.followup_evidence}) |",
        f"| contig naming | {'mismatch, harmonized' if c.chr_naming_mismatch else 'consistent'} |",
        f"| contig length conflicts | {'; '.join(c.length_conflicts) or 'none'} |",
        f"| shared contigs | {_fmt(len(c.shared_contigs))} |",
        f"| sample identity | {c.sample_note or 'not comparable'} |",
        f"| comparison permitted | {c.compatible} |",
        "",
    ]
    if c.warnings:
        out += ["**Warnings**", ""] + [f"- {w}" for w in c.warnings] + [""]
    out += [
        "## 3. Comparison settings", "",
        "| setting | value |", "| --- | --- |",
        f"| comparison key | normalized (contig, pos, ref, alt); genotype excluded |",
        f"| normalization level | `{res.normalization_level}` |",
    ]
    for k, v in res.options.items():
        out.append(f"| {k} | {v} |")
    out += [
        "",
        "## 4. Results", "",
        "| category | variants |", "| --- | ---: |",
        f"| reported in both files | {_fmt(res.n_shared)} |",
        f"| reported in the baseline file but not reported in the follow-up file | {_fmt(res.n_baseline_only)} |",
        f"| reported in the follow-up file but not reported in the baseline file | {_fmt(res.n_followup_only)} |",
        f"| reported in both files with a differing genotype | {_fmt(res.n_genotype_differences)} |",
        "",
        f"Total alleles reported in the baseline file: {_fmt(res.n_reported_baseline)}; "
        f"in the follow-up file: {_fmt(res.n_reported_followup)}.",
        "",
        "### 4.1 Reported in the baseline file but not in the follow-up file", "",
        _context_table(res.baseline_only_contexts, "follow-up"), "",
        "### 4.2 Reported in the follow-up file but not in the baseline file", "",
        _context_table(res.followup_only_contexts, "baseline"), "",
    ]
    if res.regions_used:
        out += _regions_section(res)
    if res.reference_path:
        out += _reference_section(res)
    out += [
        "## 5. How to read these numbers", "",
        CAVEAT_TEXT.replace("\n", " "), "",
        (LEFT_ALIGNED_CAVEAT if res.reference_path else NORMALIZATION_CAVEAT).replace("\n", " "),
        "",
    ]
    if res.regions_used:
        out += [REGION_CAVEAT.replace("\n", " "), ""]
    if res.contigs_baseline_only or res.contigs_followup_only:
        out += [
            "**Contig coverage.** Contigs present in only one input were not compared "
            "against anything: "
            f"baseline-only {res.contigs_baseline_only or 'none'}, "
            f"follow-up-only {res.contigs_followup_only or 'none'}. Variants on them are "
            "listed as unique to their file purely because the other file does not describe "
            "that contig.", "",
        ]
    if res.output_files or artifact_links:
        out += ["## 6. Output files", ""]
        for name, path in (artifact_links or res.output_files).items():
            out.append(f"- {name}: {path}")
        out.append("")
    return "\n".join(out)


def check_wording(text: str):
    """Return any forbidden interpretive terms found in ``text`` (QA helper)."""
    low = text.lower()
    return [t for t in FORBIDDEN_INTERPRETIVE_TERMS if t in low]


def write_report(res: ComparisonResult, outdir: str, title: str = "GTT comparison report",
                 artifact_links: Optional[Dict[str, str]] = None) -> Dict[str, str]:
    """Write ``comparison_report.md`` and ``comparison_summary.json`` into ``outdir``."""
    os.makedirs(outdir, exist_ok=True)
    md_path = os.path.join(outdir, "comparison_report.md")
    with open(md_path, "w", newline="\n", encoding="utf-8") as fh:
        fh.write(render_markdown(res, title=title, artifact_links=artifact_links))
    summary = {
        "gtt_version": __version__,
        "baseline": res.baseline_header.summary(),
        "followup": res.followup_header.summary(),
        "compatibility": {k: v for k, v in res.compatibility.as_dict().items()},
        "options": res.options,
        "baseline_stats": res.baseline_stats.as_dict(),
        "followup_stats": res.followup_stats.as_dict(),
        "counts": res.summary(),
        "baseline_only_contexts": res.baseline_only_contexts,
        "followup_only_contexts": res.followup_only_contexts,
        "contigs_compared": res.contigs_compared,
        "contigs_baseline_only": res.contigs_baseline_only,
        "contigs_followup_only": res.contigs_followup_only,
        "output_files": res.output_files,
    }
    if res.reference_path:
        summary["reference"] = dict(res.reference_summary, path=res.reference_path)
        summary["ref_mismatches"] = res.ref_mismatch_rows[:100]
    if res.regions_used:
        summary["regions"] = res.region_provenance
        summary["region_mode"] = res.region_mode
        summary["not_interrogated"] = {
            "baseline": res.n_not_interrogated_baseline,
            "followup": res.n_not_interrogated_followup,
        }
    json_path = os.path.join(outdir, "comparison_summary.json")
    with open(json_path, "w", newline="\n", encoding="utf-8") as fh:
        json.dump(summary, fh, indent=2, default=str)
    return {"report_md": md_path, "summary_json": json_path}
