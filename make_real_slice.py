"""Create a head slice of the public GIAB HG002 benchmark VCF for parser tests.

The slice is a verbatim copy of the header plus the first N data records of the
public file - real public reference data, not synthetic. Run from the repo root:

    python scripts/make_real_slice.py <hg002.vcf.gz> tests/fixtures/hg002_chr1_head.vcf 50000
"""
import gzip
import sys


def main(src, dst, n_records=50000):
    n = 0
    with gzip.open(src, "rt") as fin, open(dst, "w", newline="\n", encoding="utf-8") as fout:
        for line in fin:
            if line.startswith("#"):
                fout.write(line)
                continue
            fout.write(line)
            n += 1
            if n >= int(n_records):
                break
    print(f"wrote {n} records to {dst}")


if __name__ == "__main__":
    main(*sys.argv[1:])
