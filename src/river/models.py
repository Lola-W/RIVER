"""Input records and preset definitions."""
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

PRESETS = {
    "IlluminaSeq": "bulk",
    "ppmSeq": "duplex",
    "UDSeq": "duplex",
    "NanoSeq": "duplex",
    "HiDEF-seq": "duplex",
    "PTA-Seq": "single_cell",
}


@dataclass(frozen=True, order=True)
class Allele:
    chrom: str
    pos: int
    ref: str
    alt: str


@dataclass(frozen=True)
class Origin:
    allele: Allele
    source_sample_id: str
    callset: str = ""


@dataclass(frozen=True)
class Sample:
    sample_id: str
    method: str
    type: str
    bam: Optional[Path]
    preset: str = ""
    bam_state: str = "raw"
    strand_source: str = "sam_flags"
    fwd_pattern: str = "forward"
    rev_pattern: str = "reverse"
    family_tag: str = ""

    @property
    def orientation_labels(self):
        return self.fwd_pattern, self.rev_pattern
