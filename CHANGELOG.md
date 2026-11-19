# Changelog

All notable changes to GTT. Dates are ISO (YYYY-MM-DD).

## 0.2.0 — in development

### Added (Phase 1: region filtering and reference-based left-alignment)

- `gtt/fasta.py` — indexed random access to a reference FASTA, standard library
  only. Reads a samtools `.fai` index and builds one when it is missing
  (falling back to an in-memory index if the directory is not writable).
  Handles LF and CRLF line endings, soft-masked (lowercase) sequence, and
  contig-name differences such as `chr20` vs `20` or `chrM` vs `MT`.
  Gzip-compressed FASTA is refused with instructions rather than mis-read.
- `gtt/regions.py` — BED parsing with interval merging and intersection, and
  one single place where the BED (0-based, half-open) to VCF (1-based,
  inclusive) coordinate conversion happens.
- `normalize.left_align_allele()` — full left-alignment by the algorithm of
  Tan, Abecasis & Kang (2015), *Bioinformatics* 31(13):2202-4, the algorithm
  used by `vt normalize` and `bcftools norm`. Each allele is checked against the
  reference first; a REF that disagrees is left as written, counted, and listed.
- `gtt compare --reference FASTA` — raises `normalization_level` from
  `trim_only` to `left_aligned`, writes `ref_mismatches.tsv`, and reports per
  allele how many alleles were already aligned, moved, truncated at a contig
  start, or could not be checked.
- `gtt compare --regions BED`, `--regions-baseline BED`, `--regions-followup BED`
  — restrict the comparison to the intersection of the supplied regions.
  Alleles outside it are reported as **not interrogated** and are excluded from
  the difference counts. `--emit-not-interrogated` writes them to
  `not_interrogated.tsv`.
- `gtt compare --region-mode {start,span}` — whether the normalized POS alone
  must be inside a region (default, the `bcftools -T` convention) or the whole
  REF allele.
- `gtt faidx` and `gtt regions` subcommands for inspecting a reference index and
  BED files.
- `scripts/crosscheck_bcftools.py` — optional validation of GTT's normalization
  against `bcftools norm`. bcftools is never a runtime dependency.
- Fixtures: `synthetic_reference.fa` (3 hand-written contigs with a homopolymer,
  a CA repeat and a contig-start case), `synthetic_repeat_baseline.vcf` /
  `synthetic_repeat_followup.vcf`, and five BED files. Expected results are
  hand-computed and tabulated in `tests/fixtures/README.md`.
- 68 new tests (146 total, from 78).

### Unchanged on purpose

- The core `gtt` package still has no runtime dependencies.
- Without `--reference` and `--regions`, the output is v0.1.0's: the difference
  tables are byte-identical, the JSON summary carries exactly the v0.1.0 keys,
  and no extra files are written. `tests/test_phase1_compare.py::TestV010Parity`
  enforces this.
- The conservative report vocabulary is untouched. Alleles outside the
  interrogated regions are described as *not interrogated*, which is a statement
  about the scope of the comparison and not about the variant.

## 0.1.0

- First release: streaming VCF reader with multiallelic decomposition and
  genotype projection, parsimonious trimming, reference-build fingerprinting
  with refusal on mismatch, contig-blocked comparison engine with a recorded
  reason for every difference, and Markdown/JSON reporting with an enforced
  conservative vocabulary. 78 tests.
