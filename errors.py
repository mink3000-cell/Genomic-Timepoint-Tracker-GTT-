"""Exception hierarchy for GTT."""


class GttError(Exception):
    """Base class for all GTT errors."""


class VcfFormatError(GttError):
    """The input is not a structurally valid VCF, or violates a GTT requirement."""

    def __init__(self, message, path=None, line_no=None):
        self.path = path
        self.line_no = line_no
        loc = ""
        if path is not None:
            loc = f" [{path}"
            loc += f":{line_no}]" if line_no is not None else "]"
        super().__init__(f"{message}{loc}")


class BuildMismatchError(GttError):
    """The two inputs are not on a comparable reference genome assembly."""


class UnsortedInputError(VcfFormatError):
    """Records for one contig are interleaved with another contig."""


class FastaError(GttError):
    """The reference FASTA or its .fai index cannot be used for random access."""


class BedFormatError(GttError):
    """A BED region file is not parseable."""
