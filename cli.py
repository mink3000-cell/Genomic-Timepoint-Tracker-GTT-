"""Command-line interface: ``gtt inspect`` and ``gtt compare``."""
from __future__ import annotations

import argparse
import json
import sys

from . import __version__
from .compare import compare_files
from .errors import GttError
from .fasta import IndexedFasta
from .regions import RegionSet
from .report import write_report
from .validate import check_compatibility, infer_build, validate_header
from .vcfio import VcfReader, read_header


def _cmd_inspect(args) -> int:
    hdr = read_header(args.vcf)
    build, evidence = infer_build(hdr)
    info = hdr.summary()
    info["inferred_build"] = build
    info["build_evidence"] = evidence
    info["warnings"] = validate_header(hdr)
    if args.count:
        reader = VcfReader(args.vcf, sample=args.sample, skip_star=not args.keep_star)
        for _ in reader:
            pass
        info["stats"] = reader.stats.as_dict()
    print(json.dumps(info, indent=2))
    return 0


def _cmd_regions(args) -> int:
    sets = [RegionSet.from_bed(p) for p in args.bed]
    info = {"files": [rs.summary() for rs in sets]}
    if len(sets) > 1:
        combined = sets[0]
        for rs in sets[1:]:
            combined = combined.intersect(rs)
        info["intersection"] = combined.summary()
    print(json.dumps(info, indent=2))
    return 0


def _cmd_faidx(args) -> int:
    fa = IndexedFasta(args.fasta)
    try:
        if args.fetch:
            contig, _, span = args.fetch.partition(":")
            start_s, _, end_s = span.partition("-")
            start, end = int(start_s), int(end_s or start_s)
            print(f">{contig}:{start}-{end}")
            print(fa.fetch(contig, start - 1, end))
        else:
            print(json.dumps({**fa.summary(), "contigs": fa.contigs}, indent=2))
    finally:
        fa.close()
    return 0


def _cmd_compare(args) -> int:
    res = compare_files(
        args.baseline, args.followup, outdir=args.outdir,
        sample_baseline=args.sample_baseline, sample_followup=args.sample_followup,
        chr_naming=args.chr_naming, pass_only=args.pass_only,
        require_carried=args.require_carried, skip_star=not args.keep_star,
        skip_symbolic=args.skip_symbolic, emit_shared=args.emit_shared,
        allow_build_mismatch=args.allow_build_mismatch,
        require_same_sample=args.require_same_sample,
        reference=args.reference, regions=args.regions,
        regions_baseline=args.regions_baseline, regions_followup=args.regions_followup,
        region_mode=args.region_mode, emit_not_interrogated=args.emit_not_interrogated,
    )
    paths = write_report(res, args.outdir, title=args.title)
    res.output_files.update(paths)
    print(json.dumps(res.summary(), indent=2))
    print(f"report: {paths['report_md']}", file=sys.stderr)
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="gtt",
        description="GTT - compare a baseline and a follow-up VCF and report the differences "
                    "between the two files. GTT does not infer biological change.",
    )
    p.add_argument("--version", action="version", version=f"gtt {__version__}")
    sub = p.add_subparsers(dest="command", required=True)

    i = sub.add_parser("inspect", help="parse and validate one VCF header (optionally all records)")
    i.add_argument("vcf")
    i.add_argument("--count", action="store_true", help="stream all records and report counters")
    i.add_argument("--sample", default=None)
    i.add_argument("--keep-star", action="store_true", help="keep '*' upstream-deletion alleles")
    i.set_defaults(func=_cmd_inspect)

    c = sub.add_parser("compare", help="compare a baseline VCF against a follow-up VCF")
    c.add_argument("baseline")
    c.add_argument("followup")
    c.add_argument("-o", "--outdir", required=True)
    c.add_argument("--title", default="GTT comparison report")
    c.add_argument("--sample-baseline", default=None)
    c.add_argument("--sample-followup", default=None)
    c.add_argument("--chr-naming", choices=["strict", "harmonize"], default="strict")
    c.add_argument("--pass-only", action="store_true",
                   help="count only FILTER=PASS records as reported")
    c.add_argument("--require-carried", action="store_true",
                   help="count an allele as reported only when the genotype carries it")
    c.add_argument("--keep-star", action="store_true")
    c.add_argument("--skip-symbolic", action="store_true")
    c.add_argument("--emit-shared", action="store_true",
                   help="also write the shared-variant table (can be very large)")
    c.add_argument("--allow-build-mismatch", action="store_true",
                   help="override the refusal to compare across reference builds (not advised)")
    c.add_argument("--require-same-sample", action="store_true")
    c.add_argument("--reference", default=None, metavar="FASTA",
                   help="reference FASTA for full left-alignment of indels; the .fai index "
                        "is created if missing. Without it, alleles are trimmed only.")
    c.add_argument("--regions", default=None, metavar="BED",
                   help="compare only inside these regions; variants outside are reported "
                        "as not interrogated, never as absent")
    c.add_argument("--regions-baseline", default=None, metavar="BED",
                   help="regions interrogated by the baseline file")
    c.add_argument("--regions-followup", default=None, metavar="BED",
                   help="regions interrogated by the follow-up file; with "
                        "--regions-baseline the comparison uses the intersection")
    c.add_argument("--region-mode", choices=["start", "span"], default="start",
                   help="'start': the normalized POS must be inside a region (default); "
                        "'span': every base of REF must be")
    c.add_argument("--emit-not-interrogated", action="store_true",
                   help="also write the not_interrogated.tsv table (can be very large)")
    c.set_defaults(func=_cmd_compare)

    n = sub.add_parser("regions", help="summarize a BED file, or intersect two of them")
    n.add_argument("bed", nargs="+", help="one or two BED files")
    n.set_defaults(func=_cmd_regions)

    f = sub.add_parser("faidx", help="create or read the .fai index of a reference FASTA")
    f.add_argument("fasta")
    f.add_argument("--fetch", default=None, metavar="CONTIG:START-END",
                   help="print a 1-based inclusive slice of the sequence")
    f.set_defaults(func=_cmd_faidx)
    return p


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    try:
        return args.func(args)
    except GttError as exc:
        print(f"gtt: {type(exc).__name__}: {exc}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
