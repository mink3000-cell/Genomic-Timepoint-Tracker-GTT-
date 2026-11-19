"""Cross-check GTT's left-alignment against `bcftools norm`.

bcftools is a *validation* tool here, never a runtime dependency: GTT itself
runs on the Python standard library alone. This script is how we check that
GTT's normalization agrees with the reference implementation of the same
algorithm on real data.

What it does
------------
1. runs ``bcftools norm -m -any -f <reference> <vcf>``, which splits multiallelic
   records and left-aligns every allele;
2. reads that output with GTT's own parser but with normalization switched off,
   giving the set of keys *bcftools* produced;
3. reads the original VCF with GTT's left-alignment switched on, giving the set
   of keys *GTT* produced;
4. compares the two sets and prints the disagreements.

Usage (PowerShell)::

    $env:PYTHONPATH="."
    python scripts\\crosscheck_bcftools.py --vcf slice.vcf.gz --reference GRCh38.fa

Exit codes: 0 = the two agree, 1 = disagreements found, 3 = bcftools not found.
"""
from __future__ import annotations

import argparse
import os
import shutil
import subprocess
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from gtt.fasta import IndexedFasta            # noqa: E402
from gtt.vcfio import VcfReader               # noqa: E402


def have_bcftools(exe: str = "bcftools") -> bool:
    return shutil.which(exe) is not None


def run_bcftools_norm(vcf: str, reference: str, out_vcf: str, exe: str = "bcftools") -> str:
    """Run ``bcftools norm`` and return its stderr (which carries its own counters)."""
    # -c w: warn about REF/reference disagreements instead of exiting, so that a
    # file with a mismatch can still be compared allele by allele. GTT's own
    # behaviour at such a position is to leave the allele as written and count it.
    cmd = [exe, "norm", "-m", "-any", "-c", "w", "-f", reference,
           "-o", out_vcf, "-O", "v", vcf]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        raise RuntimeError(f"bcftools failed ({proc.returncode}): {proc.stderr.strip()}")
    return proc.stderr


def keys_from(path: str, fasta=None):
    """Set of (contig, pos, ref, alt) keys, with or without GTT left-alignment."""
    out = set()
    with VcfReader(path, fasta=fasta) as reader:
        for call in reader:
            out.add((call.contig, call.pos, call.ref, call.alt))
    return out, reader.stats


def crosscheck(vcf: str, reference: str, exe: str = "bcftools", limit: int = 20):
    """Compare GTT's keys against bcftools' keys for one VCF."""
    with tempfile.TemporaryDirectory() as td:
        normed = os.path.join(td, "bcftools_norm.vcf")
        stderr = run_bcftools_norm(vcf, reference, normed, exe=exe)
        bcf_keys, _ = keys_from(normed)             # already normalized by bcftools
        fa = IndexedFasta(reference)
        try:
            gtt_keys, stats = keys_from(vcf, fasta=fa)
        finally:
            fa.close()
    only_gtt = sorted(gtt_keys - bcf_keys)
    only_bcf = sorted(bcf_keys - gtt_keys)
    return {
        "n_gtt": len(gtt_keys),
        "n_bcftools": len(bcf_keys),
        "n_shared": len(gtt_keys & bcf_keys),
        "only_gtt": only_gtt[:limit],
        "only_bcftools": only_bcf[:limit],
        "n_only_gtt": len(only_gtt),
        "n_only_bcftools": len(only_bcf),
        "agree": not only_gtt and not only_bcf,
        "bcftools_stderr": stderr.strip(),
        "gtt_stats": stats.as_dict(),
    }


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--vcf", required=True)
    p.add_argument("--reference", required=True)
    p.add_argument("--bcftools", default="bcftools")
    args = p.parse_args(argv)

    if not have_bcftools(args.bcftools):
        print("bcftools was not found on PATH.", file=sys.stderr)
        print("This check is optional: GTT does not need bcftools to run. To install it, "
              "use WSL or conda:\n    conda install -c bioconda bcftools", file=sys.stderr)
        return 3

    res = crosscheck(args.vcf, args.reference, exe=args.bcftools)
    print(f"alleles after GTT left-alignment : {res['n_gtt']}")
    print(f"alleles after bcftools norm      : {res['n_bcftools']}")
    print(f"identical keys                   : {res['n_shared']}")
    print(f"only GTT                         : {res['n_only_gtt']}")
    print(f"only bcftools                    : {res['n_only_bcftools']}")
    if res["agree"]:
        print("\nThe two normalizations agree on every allele.")
        return 0
    print("\nFirst disagreements (contig, pos, ref, alt):")
    for k in res["only_gtt"]:
        print(f"  GTT only     : {k}")
    for k in res["only_bcftools"]:
        print(f"  bcftools only: {k}")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
