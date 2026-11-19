# Phase 1 notes — region filtering and reference-based left-alignment

GTT v0.2.0, Phase 1 of 4. Everything below was run and checked; where something
was *not* run, it says so.

---

## 1. What was built, in plain language

### The problem Phase 1 solves

Two VCF files can describe **the same DNA change in different words**. Inside a
run of identical bases there is no way to tell which copy was deleted, so one
file may say "the A at position 17 is missing" and the other "the A at position
11 is missing". Both sentences describe the same genome. GTT v0.1.0 compared
the sentences, so it counted that one change as two differences — one in each
table.

Fixing this needs the reference genome, because deciding how far left a change
can slide requires knowing what the surrounding letters actually are.

The second problem is **scope**. A variant caller only makes a statement about
the stretches of genome it could actually read. Outside those stretches the file
says nothing — and *nothing* is not the same as *no variant*. If one file was
confident over a region and the other was not, every position in the gap looks
like a difference without any difference existing.

### The two features

**`--reference FASTA` — left-alignment.** GTT now slides every insertion and
deletion as far left as the reference allows, using the algorithm published by
Tan, Abecasis & Kang in 2015, which is the algorithm `vt normalize` and
`bcftools norm` use. Equivalent spellings collapse onto one representation, so
they compare equal. The report then says `normalization_level = left_aligned`
instead of `trim_only`.

Reading the reference needed a small piece of machinery: a reference FASTA is
about 3 GB, so GTT reads an **index** (`.fai`) listing where each chromosome
starts and how many letters are on each line, then seeks straight to the bases
it needs. GTT builds that index itself if it is missing, in the format samtools
uses.

Two things are reported rather than hidden:

- an indel in a run that reaches the very start of a chromosome cannot slide any
  further; it is reported at base 1 and counted as *truncated at the contig
  start*;
- if the REF column of the VCF does not match the reference sequence at that
  position, the two files are not describing the same genome there. GTT does
  **not** silently "fix" the allele. It leaves it exactly as written, counts it,
  and lists it in `ref_mismatches.tsv`.

**`--regions BED` — interrogated regions.** A BED file lists stretches of a
chromosome. GTT compares only inside them. With `--regions-baseline` and
`--regions-followup` you give one file each, and GTT compares inside the
**intersection** — a position has to be one that *both* files make a statement
about. Anything outside is reported as **not interrogated** and is kept out of
the difference counts entirely, because the question was never asked there.

### Why that wording matters

"Not interrogated" is a statement about this comparison. "Not reported" is a
statement about a file. Neither is a statement about the person's DNA. The
conservative vocabulary from v0.1.0 is unchanged, and the new tables carry the
same discipline: every row of `not_interrogated.tsv` ends with the sentence
"This is not a statement that the variant is present or absent in the other
file."

---

## 2. Commands to test it yourself

Open the project folder in VS Code, open the PowerShell terminal, and paste
these one at a time.

### Run the whole test suite

```powershell
$env:PYTHONPATH="."
python -m unittest discover -s tests -t .
```

Expected, in about a second:

```
Ran 148 tests in 0.2s

OK (skipped=14)
```

The 14 skips are expected: 12 real-data tests that need the HG002 slice
(`scripts/make_real_slice.py`) and 2 cross-check tests that need `bcftools`,
which is optional.

### See the problem, then see it fixed

```powershell
python -m gtt.cli compare tests\fixtures\synthetic_repeat_baseline.vcf tests\fixtures\synthetic_repeat_followup.vcf -o out_trim
```

Expected: `"reported_in_both": 1, "baseline_only": 5, "followup_only": 4` and
`"normalization_level": "trim_only"`. The two files describe the same three
indels, written at different anchors, and GTT counts them as six separate
differences.

```powershell
python -m gtt.cli compare tests\fixtures\synthetic_repeat_baseline.vcf tests\fixtures\synthetic_repeat_followup.vcf -o out_ref --reference tests\fixtures\synthetic_reference.fa
```

Expected: `"reported_in_both": 4, "baseline_only": 2, "followup_only": 1`,
`"normalization_level": "left_aligned"`, `"ref_mismatches": 1`. The three
indels are now recognised as shared. What is left is the two genuinely
one-sided SNVs and the one allele whose REF disagrees with the reference.

Open `out_ref\comparison_report.md` and read section 4.4, then
`out_ref\ref_mismatches.tsv`.

### Restrict to regions

```powershell
python -m gtt.cli compare tests\fixtures\synthetic_repeat_baseline.vcf tests\fixtures\synthetic_repeat_followup.vcf -o out_bed --reference tests\fixtures\synthetic_reference.fa --regions tests\fixtures\synthetic_regions.bed --emit-not-interrogated
```

Expected: `"reported_in_both": 3`, `"baseline_only": 1`, `"followup_only": 0`,
`"not_interrogated_baseline": 2`, `"not_interrogated_followup": 2`. Open
`out_bed\not_interrogated.tsv`: the SNVs at chrT:35 and chrT:37 are there, and
they are no longer in `baseline_only.tsv` / `followup_only.tsv`.

### Look at the reference index and the regions

```powershell
python -m gtt.cli faidx tests\fixtures\synthetic_reference.fa
python -m gtt.cli faidx tests\fixtures\synthetic_reference.fa --fetch chrT:9-22
python -m gtt.cli regions tests\fixtures\synthetic_regions.bed tests\fixtures\synthetic_regions_narrow.bed
```

The `--fetch` command should print `ACAAAAAAAAGTCA`: the C anchor at position
10, the eight-A homopolymer, then `GT` and the start of the CA repeat. That is
the sequence every hand-computed test expectation was derived from.

### Check that nothing from v0.1.0 changed

```powershell
python -m gtt.cli compare tests\fixtures\synthetic_baseline.vcf tests\fixtures\synthetic_followup.vcf -o out_old
```

Expected: `"reported_in_both": 7, "baseline_only": 4, "followup_only": 3` —
exactly the v0.1.0 numbers, with `"normalization_level": "trim_only"` and no new
keys in the JSON.

### Optional: cross-check against bcftools

Only if you have bcftools (WSL, or `conda install -c bioconda bcftools`):

```powershell
python scripts\crosscheck_bcftools.py --vcf tests\fixtures\synthetic_repeat_followup.vcf --reference tests\fixtures\synthetic_reference.fa
```

Expected: `The two normalizations agree on every allele.` If bcftools is not
installed the script exits with code 3 and says so, which is not a failure.

---

## 3. Development log

### Completed

| item | detail |
| --- | --- |
| `gtt/fasta.py` | `.fai` building/parsing/writing, byte-offset random access, LF and CRLF, soft-masked sequence uppercased, `chr20`/`20`/`chrM`/`MT` name aliasing, gzip refused with instructions, in-memory fallback when the directory is read-only |
| `gtt/regions.py` | BED parsing (comments, `track`/`browser`, tab or whitespace, gzip), interval merging, intersection, `contains` / `contains_span`, one single coordinate conversion |
| `normalize.left_align_allele()` | Tan/Abecasis/Kang 2015 left-shift with a REF check first; statuses `unchanged`, `shifted`, `ref_mismatch`, `truncated_at_contig_start`, `shift_limit_reached`, `contig_not_in_reference`, `not_applicable` |
| `gtt compare` | `--reference`, `--regions`, `--regions-baseline`, `--regions-followup`, `--region-mode`, `--emit-not-interrogated` |
| new subcommands | `gtt faidx`, `gtt regions` |
| new outputs | `ref_mismatches.tsv`, `not_interrogated.tsv`, report sections 4.3 and 4.4, new JSON blocks |
| tests | 78 → 148 (14 skipped: 12 real-data, 2 bcftools) |
| fixtures | `synthetic_reference.fa` (+ `.fai`), two repeat VCFs, five BED files, all documented with hand-computed expectations in `tests/fixtures/README.md` |
| docs | README sections, `CHANGELOG.md`, this file |

Design decisions worth being able to defend in an interview:

1. **Left-alignment happens before the region test**, so a variant is judged at
   its normalized position. Sliding an indel can move it across a region edge;
   doing it in the other order would make the result depend on how the file
   happened to spell the variant.
2. **`--region-mode start` is the default** — the normalized POS must be inside
   a region, which is the convention `bcftools -T` uses. `span` (the whole REF
   allele inside one interval) is stricter and available, but making it the
   default would make GTT disagree with the standard tools at region edges.
3. **A REF mismatch is never corrected.** Silently rewriting the allele to match
   the FASTA would make a wrong-build comparison look clean.
4. **Multiple BED files are intersected, not unioned.** A union would compare
   positions that only one file makes a statement about, which is exactly the
   artefact the feature exists to remove.
5. **The index is written next to the FASTA** when possible, so a second run is
   fast and samtools can use the same file.

### Known issues and limitations

1. **The bcftools cross-check has not been executed.** bcftools is not installed
   in the environment this was developed in, so `tests/test_bcftools_crosscheck.py`
   skipped and the script has never run against real bcftools output. Agreement
   between GTT and bcftools is therefore **an expectation, not a verified fact**,
   and the script itself is untested code. Running it is the first item of
   Phase 2.
2. **Left-alignment is verified against hand-computed cases only** — 14 of them,
   on a 3-contig synthetic reference. Every expectation was worked out on paper
   from the sequence before the code was run, which is a real check, but it is
   not a check against an independent implementation on real data.
3. **bgzip-compressed FASTA is not supported.** Random access into a bgzip file
   needs a `.gzi` index and a block-level decompressor. GTT refuses with an
   instruction to decompress rather than reading it wrongly. The GRCh38
   reference must be uncompressed on disk (about 3 GB).
4. **A stale `.fai` is only partly detectable.** If the index claims more
   sequence than the file holds, the short read is caught. If it is merely wrong
   in the middle of the file, GTT can only warn that the index is older than the
   FASTA.
5. **No performance measurement on a whole chromosome yet.** Left-alignment adds
   file seeks per indel; the 64-base window buffer keeps that to roughly one read
   per indel, but the real cost on chr20 of HG002 is unmeasured. Phase 2 will
   produce that number.
6. **`requires-python` was corrected to `>=3.10`.** The package already used
   `@dataclass(slots=True)`, which is 3.10+, so the previous `>=3.9` was wrong.
7. **`__version__` is now 0.2.0**, so reports generated now say v0.2.0 even when
   run with v0.1.0 behaviour. The *results* are unchanged; only the version
   string and the timestamp differ.

### Next step

Phase 2, the noise-floor experiment. The first three actions are: run the
bcftools cross-check on a real chr20 slice (issue 1 above); identify and verify
the two HG002 callsets against the official GIAB/NIST source, recording exact
URL, version and checksum before downloading anything; and write
`docs/predictions.md` before looking at any result.

**Environment note.** The machine this was developed on currently has no network
access — a DNS lookup for the NCBI FTP host failed — so no GIAB file has been
downloaded or verified. If that is still true at the start of Phase 2, you will
be given the exact PowerShell commands to download and checksum the files on
your own machine, and the analysis will run against the files you provide.
