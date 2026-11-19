# GTT — Genomic Timepoint Tracker

GTT compares two processed VCF files — a **baseline** and a **follow-up**
observation — and produces a transparent difference report: variants reported in
both files, variants reported only in the baseline file, and variants reported
only in the follow-up file.

**What GTT does not do.** It does not infer biological change. A variant absent
from a VCF is reported as *not reported in that file*, never as a confirmed
loss; a variant present only in the follow-up file is reported as *reported only
in the follow-up file*, never as a newly acquired or treatment-induced mutation.
Those words are blocked at the report layer, not left to the writer's
discretion. External annotation (VEP, ClinVar, population frequencies, Open
Targets, ACMG-style evidence) is a possible future extension and is **not** a
prerequisite for the core workflow.

Zero runtime dependencies — Python standard library only.

## Install

```
pip install -e .
```

Or run from a checkout without installing:

```
PYTHONPATH=. python -m gtt.cli --help      # Windows: $env:PYTHONPATH="."
```

## Use

```
gtt inspect sample.vcf.gz --count
gtt compare baseline.vcf.gz followup.vcf.gz -o results/
```

`compare` writes, into `-o`:

| file | contents |
| --- | --- |
| `comparison_report.md` | the human-readable report, with inputs, build checks, counts and caveats |
| `comparison_summary.json` | the same numbers in machine-readable form |
| `baseline_only.tsv` | variants reported in the baseline file but not in the follow-up file |
| `followup_only.tsv` | variants reported in the follow-up file but not in the baseline file |
| `genotype_differences.tsv` | variants reported in both files whose genotype differs |
| `shared.tsv` | all shared variants (only with `--emit-shared`; can be millions of rows) |

Useful options: `--pass-only` (count only `FILTER=PASS` records as reported),
`--require-carried` (require the genotype to carry the allele),
`--chr-naming harmonize` (compare `chr1` against `1`), `--sample-baseline` /
`--sample-followup` (pick a sample in a multi-sample VCF), `--keep-star`,
`--skip-symbolic`.

### Reference-based left-alignment (v0.2.0)

```
gtt compare baseline.vcf.gz followup.vcf.gz -o results/ --reference GRCh38.fa
```

Without `--reference`, alleles are only trimmed, so the same indel written at
two different anchors inside a repeat appears once in each difference table.
With `--reference`, every allele is fully left-aligned by the algorithm of Tan,
Abecasis & Kang (2015) — the algorithm behind `vt normalize` and `bcftools
norm` — and the two spellings collapse onto one key. The report then states
`normalization_level = left_aligned`.

The FASTA must be uncompressed. Its `.fai` index is created if missing (GTT
writes a samtools-compatible index; if the directory is read-only the index is
kept in memory for that run). Every REF allele is checked against the reference
sequence: a REF that disagrees is **never corrected**. It is left exactly as the
file wrote it, counted in the report, and listed in `ref_mismatches.tsv`,
because a disagreement means the VCF and the FASTA are not describing the same
genome at that position.

```
gtt faidx GRCh38.fa                       # build/read the index, print contigs
gtt faidx GRCh38.fa --fetch chr20:100-130 # print a slice, 1-based inclusive
```

### Region restriction

```
gtt compare baseline.vcf.gz followup.vcf.gz -o results/ `
    --regions-baseline baseline_confident.bed `
    --regions-followup followup_confident.bed
```

A variant caller only makes a statement about the regions it could examine.
`--regions BED` restricts the comparison to those regions; `--regions-baseline`
and `--regions-followup` describe one input each, and when more than one BED is
given the comparison region is their **intersection** — a position must be
interrogated in every supplied file to be compared.

Alleles outside the comparison region are reported as **not interrogated**.
They are not counted as baseline-only or follow-up-only, because the question
was never asked there. `--emit-not-interrogated` writes them to
`not_interrogated.tsv`. `--region-mode start` (default) tests the normalized
POS, the convention `bcftools -T` uses; `--region-mode span` requires every base
of REF to be inside one interval.

```
gtt regions a.bed b.bed    # interval/base counts for each file and their intersection
```

Left-alignment is applied **before** the region test, so a variant is judged at
its normalized position.

## Module structure

| module | responsibility |
| --- | --- |
| `gtt/model.py` | `AlleleCall` — one row per (record, ALT allele); comparison and site keys; genotype-status vocabulary |
| `gtt/normalize.py` | parsimonious allele trimming, reference-based left-alignment, variant typing |
| `gtt/fasta.py` | indexed random access to a reference FASTA; builds a samtools-compatible `.fai` |
| `gtt/regions.py` | BED parsing, interval merging, intersection, and the BED/VCF coordinate conversion |
| `gtt/vcfio.py` | streaming reader: header metadata, record parsing, multiallelic decomposition, genotype projection, structural validation |
| `gtt/validate.py` | reference-build fingerprinting, contig/naming/sample compatibility, refusal logic |
| `gtt/compare.py` | contig-blocked comparison engine producing the three variant sets and their explanations |
| `gtt/report.py` | Markdown + JSON rendering with the enforced conservative vocabulary |
| `gtt/cli.py` | `gtt inspect`, `gtt compare`, `gtt faidx`, `gtt regions` |

## Design decisions

**Reference genome compatibility.** The assembly is fingerprinted from contig
lengths (chr1 = 248,956,422 → GRCh38; 249,250,621 → GRCh37; 248,387,328 →
T2T-CHM13v2.0), cross-checked against `##reference` and `assembly=`. A build
mismatch is **refused**, not warned about: GTT has no liftover, and comparing
GRCh37 coordinates against GRCh38 coordinates produces a report that is almost
entirely spurious differences. `--allow-build-mismatch` exists as an explicit,
logged override. A `chr1`-vs-`1` naming mismatch is a separate, recoverable
condition (`--chr-naming harmonize`). Conflicting lengths for a shared contig
are a hard error.

**Variant normalization.** Each ALT is normalized independently after
decomposition: uppercase, trim the shared suffix (leaving ≥ 1 base), then the
shared prefix (leaving ≥ 1 base, advancing POS). This matches the vt/bcftools
right-then-left convention, so `CTTA>CA` and `CTT>C` compare equal.
Symbolic alleles, breakends and `*` are never trimmed. Without a reference this
is **not** full left-alignment, and the report says so
(`normalization_level = trim_only`): the same indel anchored at different
positions inside a homopolymer appears once in each difference table instead of
as a shared variant. Supplying `--reference` adds the left-shifting step and
raises the level to `left_aligned`. Two limits remain and are both counted in
the report: an indel whose leftmost representation would start before base 1 is
reported at base 1 (`truncated_at_contig_start`), and an allele whose REF
disagrees with the reference is left as written (`ref_mismatch`).

**Multiallelic records.** Split to one row per ALT, each retaining
`multiallelic`, `n_alts` and its original ALT index, so a two-ALT site that
matches on one allele reports one shared and one unique allele rather than a
whole-site mismatch. `Number=A`/`R` semantics are preserved by projecting the
genotype onto the specific ALT index.

**Genotype fields.** `FORMAT` is indexed by key; `GT` is parsed for phased,
unphased, haploid and missing calls. Genotype is **not** part of the comparison
key. A variant reported het in the baseline and hom in the follow-up is a
*shared* variant with a reported genotype difference — a different statement
from presence or absence. Agreement is judged per allele (zygosity status and
copy number), so `1/2` and `0/1` can both report one copy of the same ALT.

**Missing or unreported variants.** Every unique variant carries the reason the
other file does not report it: `no_record_at_site`, `site_record_other_allele`,
`allele_present_but_filtered`, `allele_present_but_no_call`, or
`allele_present_but_hom_ref`. For a benchmark VCF, absence outside the
high-confidence region means *not interrogated*, which the report says rather
than implying absence of the variant. With `--regions` that distinction becomes
structural: such alleles are removed from the comparison entirely and counted
separately, so they can never be mistaken for a difference.

**Memory and ordering.** The engine holds one contig block from each file at a
time, so peak memory scales with the largest chromosome rather than the genome.
The two files may list contigs in different orders; a contig that reappears
after another has begun raises `UnsortedInputError` rather than silently
producing a wrong answer.

## Tests

```
PYTHONPATH=. python -m unittest discover -s tests -t .
```

Synthetic fixtures in `tests/fixtures/` carry a `##gtt_synthetic` header, which
GTT surfaces as a warning and as a banner on any report built from them. They
are test data with hand-computed expected results — not biological observations.
`tests/fixtures/README.md` lists the expected comparison outcome variant by
variant. Real-data tests run against a verbatim head slice of the public Genome
in a Bottle HG002 GRCh38 NISTv4.2.1 benchmark VCF and skip when that slice has
not been generated (`scripts/make_real_slice.py`).

`scripts/crosscheck_bcftools.py` compares GTT's left-alignment against
`bcftools norm` on a real slice. bcftools is used for **validation only** and is
never a runtime dependency; `tests/test_bcftools_crosscheck.py` skips when it is
not installed.
