"""Derive a perturbed copy of a real VCF slice, with a manifest of expected results.

This produces TEST DATA, not an observation. The output carries a
``##gtt_synthetic`` header so GTT flags it in every report generated from it.
The perturbations are mechanical edits to a copy of a public reference file;
they do not model mutation, selection, treatment response, or any biological
process whatsoever.

Perturbations applied (deterministic, index-based):

* every ``--drop-every``-th record is removed from the follow-up copy;
* every ``--filter-every``-th record keeps its alleles but gets FILTER=LowQual;
* every ``--genotype-every``-th record with GT 0/1 is rewritten to 1/1.

The manifest records the normalized allele keys affected by each perturbation so
a test can assert GTT's output against the recipe rather than against itself.

    python scripts/make_perturbed_slice.py tests/fixtures/hg002_chr1_head.vcf \
        tests/fixtures/hg002_chr1_head_perturbed.vcf \
        tests/fixtures/hg002_chr1_head_perturbed.manifest.json
"""
from __future__ import annotations

import argparse
import json

from gtt.model import format_key
from gtt.normalize import NON_VARIANT_ALTS, is_symbolic, normalize_allele

SYNTHETIC_NOTE = (
    "mechanically perturbed copy of a public GIAB HG002 slice, created to test the GTT "
    "comparison engine; not a biological observation and not a second timepoint"
)


def _keys(chrom, pos, ref, alt_field, skip_star=True):
    out = []
    for alt in alt_field.split(","):
        alt = alt.upper()
        if alt in NON_VARIANT_ALTS and (alt != "*" or skip_star):
            continue
        if is_symbolic(alt):
            continue
        npos, nref, nalt = normalize_allele(int(pos), ref.upper(), alt)
        out.append(format_key((chrom, npos, nref, nalt)))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("src")
    ap.add_argument("dst")
    ap.add_argument("manifest")
    ap.add_argument("--drop-every", type=int, default=500)
    ap.add_argument("--filter-every", type=int, default=700)
    ap.add_argument("--genotype-every", type=int, default=900)
    args = ap.parse_args(argv)

    dropped, filtered, genotyped = [], [], []
    n_records = 0
    with open(args.src, encoding="utf-8") as fin, \
            open(args.dst, "w", newline="\n", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("##"):
                fout.write(line)
                continue
            if line.startswith("#CHROM"):
                fout.write(f"##gtt_synthetic={SYNTHETIC_NOTE}\n")
                fout.write(line)
                continue
            n_records += 1
            fields = line.rstrip("\n").split("\t")
            chrom, pos, ref, alt = fields[0], fields[1], fields[3], fields[4]
            if n_records % args.drop_every == 0:
                dropped.extend(_keys(chrom, pos, ref, alt))
                continue
            if n_records % args.filter_every == 0:
                fields[6] = "LowQual"
                filtered.extend(_keys(chrom, pos, ref, alt))
            if n_records % args.genotype_every == 0 and len(fields) > 9:
                fmt = fields[8].split(":")
                if "GT" in fmt:
                    i = fmt.index("GT")
                    vals = fields[9].split(":")
                    if i < len(vals) and vals[i] in ("0/1", "0|1", "1|0"):
                        vals[i] = "1/1"
                        fields[9] = ":".join(vals)
                        genotyped.extend(_keys(chrom, pos, ref, alt))
            fout.write("\t".join(fields) + "\n")

    manifest = {
        "source": args.src,
        "output": args.dst,
        "note": SYNTHETIC_NOTE,
        "n_source_records": n_records,
        "drop_every": args.drop_every,
        "filter_every": args.filter_every,
        "genotype_every": args.genotype_every,
        "expected_baseline_only_keys": sorted(set(dropped)),
        "expected_filtered_keys": sorted(set(filtered)),
        "expected_genotype_changed_keys": sorted(set(genotyped)),
    }
    with open(args.manifest, "w", newline="\n", encoding="utf-8") as fh:
        json.dump(manifest, fh, indent=2)
    print(json.dumps({k: v for k, v in manifest.items()
                      if not k.startswith("expected")}, indent=2))
    print(f"expected baseline-only alleles: {len(manifest['expected_baseline_only_keys'])}")


if __name__ == "__main__":
    main()
