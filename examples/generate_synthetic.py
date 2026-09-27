#!/usr/bin/env python3
"""Create a deterministic toy reference and BAMs, following tests/test_river.py.

The 500-base chr1 is artificial, not a segment of GRCh38. No study data,
downloads, random seeds, or imports from the test suite are required.
"""
import argparse
from pathlib import Path
import shutil

import pysam


def record(name, pos, alt=False, reverse=False, tags=()):
    """Make a 100-base alignment spanning a 1-based candidate position."""
    read = pysam.AlignedSegment()
    read.query_name = name
    read.flag = 16 if reverse else 0
    read.reference_id = 0
    read.reference_start = pos - 40
    read.mapping_quality = 60
    read.cigarstring = "100M"
    bases = list("A" * 100)
    if alt:
        bases[39] = "G"
    read.query_sequence = "".join(bases)
    read.query_qualities = pysam.qualitystring_to_array("I" * 100)
    read.set_tags(list(tags))
    return read


def generate(outdir):
    # A new directory prevents accidentally replacing the user's inputs.
    outdir.mkdir(parents=True, exist_ok=False)
    reference = outdir / "reference.fa"
    reference.write_text(">chr1\n" + "A" * 500 + "\n", encoding="ascii")
    pysam.faidx(str(reference))
    bamdir = outdir / "bams"
    bamdir.mkdir()
    header = {"HD": {"VN": "1.6", "SO": "coordinate"},
              "SQ": [{"SN": "chr1", "LN": 500}]}
    for sample in ("bulk", "duplex", "cell_01", "cell_02"):
        reads = []
        for pos in (50, 200):
            is_cell = sample.startswith("cell_")
            count = 1 if is_cell else 10
            alt_count = int(pos == 50) if is_cell else 2
            tags = [("st", "MIXED"), ("et", "MIXED")] if sample == "duplex" else []
            reads.extend(record(f"{sample}_{pos}_{i}", pos, alt=i < alt_count,
                                reverse=bool(i % 2), tags=tags) for i in range(count))
        if sample == "duplex":
            # Raw ppmSeq selection must discard this extra ALT observation.
            reads.append(record("duplex_unmixed", 50, alt=True,
                                tags=[("st", "MIXED"), ("et", "PLUS")]))
        path = bamdir / f"{sample}.bam"
        with pysam.AlignmentFile(str(path), "wb", header=header) as output:
            for read in sorted(reads, key=lambda r: r.reference_start):
                output.write(read)
        pysam.index(str(path))
    for name in ("samples.tsv", "variants.tsv"):
        shutil.copyfile(Path(__file__).with_name(name), outdir / name)
    print(f"Created synthetic inputs: {outdir.resolve()}")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=Path("examples/demo/data"),
                        help="New input directory (default: examples/demo/data)")
    args = parser.parse_args()
    try:
        generate(args.outdir)
    except FileExistsError:
        parser.exit(1, f"Input directory already exists; choose a new --outdir: {args.outdir}\n")


if __name__ == "__main__":
    main()
