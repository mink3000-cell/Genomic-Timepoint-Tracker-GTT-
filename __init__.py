"""GTT - Genomic Timepoint Tracker.

Compares two processed VCF files representing a baseline and a follow-up
observation of the same subject and produces a transparent difference report.

GTT reports what is *recorded in the files*. It does not infer biological
gain or loss of variants. Absence of a variant from a VCF is reported as
"not reported in <file>", never as a confirmed biological change.
"""
__version__ = "0.2.0"

from .errors import GttError, VcfFormatError, BuildMismatchError, FastaError, BedFormatError
from .model import AlleleCall, variant_key
from .vcfio import VcfHeader, VcfReader, open_vcf, read_header
from .normalize import normalize_allele, classify_variant, left_align_allele
from .fasta import IndexedFasta, build_fai, read_fai, write_fai
from .regions import RegionSet, build_comparison_regions
from .validate import infer_build, validate_header, check_compatibility, CompatibilityReport
from .compare import compare_files, ComparisonResult
from .report import write_report

__all__ = [
    "__version__", "GttError", "VcfFormatError", "BuildMismatchError",
    "FastaError", "BedFormatError",
    "AlleleCall", "variant_key", "VcfHeader", "VcfReader", "open_vcf", "read_header",
    "normalize_allele", "classify_variant", "left_align_allele",
    "IndexedFasta", "build_fai", "read_fai", "write_fai",
    "RegionSet", "build_comparison_regions",
    "infer_build", "validate_header",
    "check_compatibility", "CompatibilityReport", "compare_files",
    "ComparisonResult", "write_report",
]
