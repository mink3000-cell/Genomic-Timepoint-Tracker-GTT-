# Test fixtures

## Synthetic fixtures (NOT biological observations)

Every file in this directory whose name starts with `synthetic_` or `malformed_`
is hand-written test data created to exercise GTT's parsing, validation and
comparison logic. Each carries a `##gtt_synthetic=` header line, which GTT
detects and surfaces as a warning and as a banner on any report generated from
it. These files do not represent any individual, any sequencing run or any
observed genomic change, and no biological conclusion may be drawn from them.

`synthetic_baseline.vcf` vs `synthetic_followup.vcf` - expected results under
default options (`pass_only=False`, `require_carried=False`, `skip_star=True`):

| category | count | members |
| --- | ---: | --- |
| reported in both files | 7 | chr1-100000-A-G, chr1-200000-C-T, chr1-300000-CTT-C, chr1-800000-A-G, chr1-900000-T-C, chr2-500000-G-A, chr2-800000-G-A |
| reported in baseline only | 4 | chr1-400000-T-A, chr1-450000-G-C, chr2-500000-G-T, chrX-1000000-A-T |
| reported in follow-up only | 3 | chr1-450000-G-T, chr1-600000-C-G, chr2-700000-A-AT |
| shared with a differing genotype | 2 | chr1-200000-C-T (het vs hom), chr1-900000-T-C (het vs no-call) |

Cases covered: identical calls; same allele with a different genotype; the same
deletion written parsimoniously (`CTT>C`) in one file and non-parsimoniously
(`CTTA>CA`) in the other; a multiallelic record where one of two ALTs is shared;
two different ALT alleles at the same position; a `*` upstream-deletion
placeholder; a haploid genotype; a contig present in one file only; a record
filtered in one file (`LowQual`); and a no-call genotype.

## Reference-based normalization fixtures (v0.2.0)

`synthetic_reference.fa` is a 3-contig hand-written sequence. It is not a real
genome and no position in it corresponds to any position in any real assembly.
Its `.fai` index ships with it so the test suite does not write into the source
tree; delete the `.fai` and it is rebuilt on the next run.

```
chrT (40 bp)  ACGTACGTAC AAAAAAAA GT CACACACACACA TTGGCCAA
pos           1          11       19 21           33
                         ^ 8-base A homopolymer   ^ six CA units
              anchor base for the homopolymer: C at position 10
              anchor base for the CA repeat  : T at position 20
chrS (10 bp)  AAAAACCCCC      - an A run with no base before it, for contig-start tests
chrL (20 bp)  acgtacgtacgtACGTACGT - soft-masked, for case-insensitivity tests
```

`synthetic_repeat_baseline.vcf` vs `synthetic_repeat_followup.vcf` write the
same three indels at different anchors. Hand-computed expected results:

| allele as written | left-aligned key | note |
| --- | --- | --- |
| baseline chrT:17 `AA>A`, follow-up chrT:11 `AA>A` | `chrT-10-CA-C` | one A removed from the homopolymer |
| baseline chrT:24 `ACA>A`, follow-up chrT:30 `ACA>A` | `chrT-20-TCA-T` | one CA unit removed from the repeat |
| baseline chrS:3 `AA>A`, follow-up chrS:2 `AA>A` | `chrS-1-AA-A` | cannot shift past base 1; counted as truncated at the contig start |
| baseline chrT:5 `C>T` | not aligned | REF disagrees with the reference (`A`); counted and listed, never corrected |

| run | reported in both | baseline only | follow-up only |
| --- | ---: | ---: | ---: |
| no `--reference` (trim only) | 1 | 5 | 4 |
| `--reference synthetic_reference.fa` | 4 | 2 | 1 |
| `--reference` + `--regions synthetic_regions.bed` | 3 | 1 | 0 (plus 2 + 2 not interrogated) |

BED fixtures, all 0-based half-open as BED requires:

| file | spans | purpose |
| --- | --- | --- |
| `synthetic_regions.bed` | chrT 0-30 | the main region test; carries a comment and a `track` line |
| `synthetic_regions_narrow.bed` | chrT 0-15, chrS 0-10 | intersection test (chrT 0-30 AND chrT 0-15 = 0-15) |
| `synthetic_regions_edge.bed` | chrT 0-20 | `--region-mode span` test: `TCA` starting at POS 20 spills over the edge |
| `synthetic_regions_nochr.bed` | T 0-30 | contig naming without the `chr` prefix |
| `synthetic_regions_split.bed` | chrT 0-10, 10-20, 18-30 | touching and overlapping intervals must merge into one |

## Real public data

`hg002_chr1_head.vcf` (created by the development notebook, not checked in) is a
verbatim head slice of the Genome in a Bottle HG002 GRCh38 NISTv4.2.1 benchmark
VCF. It is real public reference data, used for parser and throughput testing.
Any comparison built from it by perturbation is labelled as derived test data,
not as a second observation of the subject.
