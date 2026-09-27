"""RIVER command-line workflow."""
import argparse
import csv
import hashlib
import json
import shutil
import subprocess
import sys
from dataclasses import asdict
from datetime import datetime, timezone
from pathlib import Path

import pysam

from . import __version__
from .inputs import load_samples, load_variants, validate_alignments
from .pileup import MPILEUP_OPTIONS, SUPPORT_THRESHOLD, WILSON_Z, query_sample
from .selection import prepare_bam
from .support import summarize

READ_FIELDS = ["chrom", "pos", "ref", "alt", "sample_id", "method", "type", "ref_count", "alt_count", "depth", "vaf", "lower_ci", "upper_ci", "qualifying_support", "query_status", "coverage_status", "ref_forward", "ref_reverse", "ref_unclassified", "alt_forward", "alt_reverse", "alt_unclassified"]
VARIANT_FIELDS = ["chrom", "pos", "ref", "alt", "source_sample_id", "source_method", "callset", "supporting_samples", "supporting_methods", "supporting_cells_by_method", "supporting_cell_counts_by_method", "supporting_cell_classes_by_method", "other_supporting_methods", "n_other_supporting_methods", "other_method_support_class"]
SUMMARY_FIELDS = ["source_sample_id", "source_method", "callset", "n_variants", "n_support_0", "n_support_1", "n_support_ge2", "proportion_support_0", "proportion_support_1", "proportion_support_ge2"]


def _utc():
    return datetime.now(timezone.utc).isoformat()


def _write_json(path, value):
    temporary = path.with_suffix(path.suffix + ".tmp")
    with temporary.open("w") as handle:
        json.dump(value, handle, indent=2, default=str, allow_nan=False)
        handle.write("\n")
    temporary.replace(path)


def _field(value):
    if value is None:
        return "NA"
    if isinstance(value, bool):
        return int(value)
    if isinstance(value, (dict, list)):
        return json.dumps(value, separators=(",", ":"), sort_keys=True)
    return value


def _writer(handle, fields):
    writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
    writer.writeheader()
    return writer


def _write_table(path, fields, rows):
    with path.open("x", newline="") as handle:
        writer = _writer(handle, fields)
        for row in rows:
            writer.writerow({key: _field(value) for key, value in row.items()})


def _metadata(path, digest=False):
    path = Path(path).resolve()
    info = {"path": str(path)}
    if path.is_file():
        stat = path.stat()
        info.update(size_bytes=stat.st_size, modified_ns=stat.st_mtime_ns)
        if digest:
            info["sha256"] = hashlib.sha256(path.read_bytes()).hexdigest()
    return info


def run(args):
    if args.threads < 1:
        raise ValueError("--threads must be at least 1")
    outdir = args.outdir.resolve()
    try:
        outdir.mkdir(parents=True, exist_ok=False)
    except FileExistsError as exc:
        raise ValueError(f"Output directory already exists; choose a new --outdir: {outdir}") from exc
    metadata_path = outdir / "run.json"
    metadata = {"status": "running", "started_at": _utc(), "river_version": __version__,
                "python": sys.version, "pysam": pysam.__version__,
                "inputs": {"variants": _metadata(args.variants, True), "samples": _metadata(args.samples, True), "reference": _metadata(args.reference)},
                "settings": {"coordinate_system": "GRCh38; 1-based variants", "threads": args.threads,
                             "mpileup_options": MPILEUP_OPTIONS, "wilson_z": WILSON_Z,
                             "bulk_duplex_lower_ci_strictly_greater_than": SUPPORT_THRESHOLD,
                             "single_cell_min_alt_count": 1, "method_combination": "any qualifying dataset",
                             "origin_exclusion": "entire method", "overlapping_mate_handling": "samtools default"},
                "samples": [], "preparation": {}, "commands": []}
    _write_json(metadata_path, metadata)
    try:
        samples = load_samples(args.samples)
        origins = load_variants(args.variants, samples)
        alleles = list(dict.fromkeys(origin.allele for origin in origins))
        reference = args.reference.resolve()
        metadata["inputs"]["reference_index"] = _metadata(Path(str(reference) + ".fai"), True)
        for sample in samples.values():
            sample_metadata = asdict(sample)
            if sample.bam:
                sample_metadata["bam_metadata"] = _metadata(sample.bam)
            metadata["samples"].append(sample_metadata)
        lengths = validate_alignments(reference, samples, alleles)
        metadata["candidate_alleles"] = len(alleles)
        metadata["unique_origins"] = len(origins)
        metadata["reference_contig_lengths"] = {chrom: lengths[chrom] for chrom in sorted({allele.chrom for allele in alleles})}
        samtools = shutil.which(args.samtools)
        if not samtools:
            raise ValueError(f"SAMtools executable not found: {args.samtools}")
        result = subprocess.run([samtools, "--version"], check=True, text=True, capture_output=True)
        metadata["samtools"] = {"path": samtools, "version": result.stdout.splitlines()[:2]}
        workdir = outdir / "intermediate"
        workdir.mkdir()
        bed = workdir / "candidate_sites.bed"
        # The BED is zero-based; only candidate positions are returned by mpileup.
        contig_order = {chrom: index for index, chrom in enumerate(lengths)}
        with bed.open("x") as handle:
            for chrom, pos in sorted({(allele.chrom, allele.pos) for allele in alleles}, key=lambda item: (contig_order[item[0]], item[1])):
                handle.write(f"{chrom}\t{pos - 1}\t{pos}\n")
        _write_json(metadata_path, metadata)
        read_path = outdir / "read_support.tsv.tmp"
        read_count = 0
        with read_path.open("x", newline="") as read_handle:
            writer = _writer(read_handle, READ_FIELDS)

            def rows():
                nonlocal read_count
                for index, sample in enumerate(samples.values()):
                    if sample.bam is None:
                        continue
                    print(f"Querying {sample.sample_id} ({sample.method})", file=sys.stderr, flush=True)
                    bam = sample.bam
                    if sample.type == "duplex" and sample.bam_state == "raw":
                        # Numeric names prevent input identifiers becoming file paths.
                        bam = workdir / f"sample_{index + 1:04d}.selected.bam"
                        stats = prepare_bam(sample, bam, samtools, args.threads)
                        metadata["preparation"][sample.sample_id] = {"output_bam": str(bam), **stats}
                        _write_json(metadata_path, metadata)
                    for row in query_sample(sample, bam, reference, alleles, bed, samtools, metadata["commands"]):
                        writer.writerow({key: _field(value) for key, value in row.items()})
                        read_count += 1
                        yield row

            variants, summaries = summarize(origins, samples, rows())
        expected = len(alleles) * sum(sample.bam is not None for sample in samples.values())
        if read_count != expected:
            raise RuntimeError(f"Incomplete query matrix: expected {expected} rows, received {read_count}")
        _write_table(outdir / "variant_support.tsv.tmp", VARIANT_FIELDS, variants)
        _write_table(outdir / "summary.tsv.tmp", SUMMARY_FIELDS, summaries)
        for name in ("read_support.tsv", "variant_support.tsv", "summary.tsv"):
            (outdir / (name + ".tmp")).replace(outdir / name)
        metadata.update(status="complete", finished_at=_utc(), read_support_rows=read_count,
                        variant_support_rows=len(variants), summary_rows=len(summaries),
                        outputs=["read_support.tsv", "variant_support.tsv", "summary.tsv"])
        _write_json(metadata_path, metadata)
        print(f"RIVER completed: {len(alleles)} alleles, {len(origins)} origins -> {outdir}")
        return 0
    except BaseException as exc:
        metadata.update(status="failed", finished_at=_utc(), error=f"{type(exc).__name__}: {exc}")
        _write_json(metadata_path, metadata)
        raise


def main(argv=None):
    parser = argparse.ArgumentParser(description="RIVER: cross-method read support for candidate SNVs")
    parser.add_argument("--version", action="version", version=__version__)
    subparsers = parser.add_subparsers(dest="command", required=True)
    command = subparsers.add_parser("run", help="Query alleles and summarize support from other sequencing methods")
    command.add_argument("--variants", required=True, type=Path)
    command.add_argument("--samples", required=True, type=Path)
    command.add_argument("--reference", required=True, type=Path)
    command.add_argument("--outdir", required=True, type=Path)
    command.add_argument("--samtools", default="samtools")
    command.add_argument("--threads", type=int, default=1)
    args = parser.parse_args(argv)
    try:
        return run(args)
    except (ValueError, OSError, RuntimeError, subprocess.SubprocessError) as exc:
        print(f"RIVER error: {exc}", file=sys.stderr)
        return 1
