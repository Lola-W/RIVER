"""Validate the two input tables and their alignment/reference dependencies."""
import csv
import re
from pathlib import Path

import pysam

from .models import Allele, Origin, PRESETS, Sample


def _rows(path, required, allowed):
    path = Path(path)
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        names = reader.fieldnames or []
        if len(names) != len(set(names)):
            raise ValueError(f"{path}: duplicate column names")
        missing = set(required) - set(names)
        unknown = set(names) - set(allowed)
        if missing:
            raise ValueError(f"{path}: missing columns: {', '.join(sorted(missing))}")
        if unknown:
            raise ValueError(f"{path}: unknown columns: {', '.join(sorted(unknown))}")
        for number, row in enumerate(reader, 2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}:{number}: inconsistent number of TSV fields")
            row = {key: value.strip() for key, value in row.items()}
            if any(not row.get(key) for key in required):
                raise ValueError(f"{path}:{number}: empty required field")
            yield number, row


def _validate_orientation(source, fwd, rev, location):
    if source == "sam_flags":
        if (fwd, rev) not in (("forward", "reverse"), ("F1R2", "F2R1")):
            raise ValueError(f"{location}: sam_flags requires forward/reverse or F1R2/F2R1")
    else:
        if source != "read_name" and not re.fullmatch(r"tag:[A-Za-z][A-Za-z0-9]", source):
            raise ValueError(f"{location}: invalid strand_source {source!r}")
        if not fwd or not rev or fwd == rev:
            raise ValueError(f"{location}: distinct forward and reverse patterns are required")
        try:
            re.compile(fwd)
            re.compile(rev)
        except re.error as exc:
            raise ValueError(f"{location}: invalid strand pattern: {exc}") from exc


def load_samples(path):
    path = Path(path).resolve()
    required = {"sample_id", "method", "type", "bam"}
    optional = {"preset", "bam_state", "strand_source", "fwd_pattern", "rev_pattern", "family_tag"}
    samples = {}
    method_types = {}
    bam_paths = set()
    for number, row in _rows(path, required, required | optional):
        location = f"{path}:{number}"
        sid, method, kind = row["sample_id"], row["method"], row["type"]
        if sid in samples:
            raise ValueError(f"{location}: duplicate sample_id {sid!r}")
        if kind not in {"bulk", "duplex", "single_cell"}:
            raise ValueError(f"{location}: type must be bulk, duplex or single_cell")
        if method in method_types and method_types[method] != kind:
            raise ValueError(f"{location}: method {method!r} has inconsistent types")
        method_types[method] = kind
        preset = row.get("preset") or (method if method in PRESETS else "")
        if method in PRESETS and PRESETS[method] != kind:
            raise ValueError(f"{location}: {method} requires type={PRESETS[method]}")
        if preset and preset not in PRESETS:
            raise ValueError(f"{location}: unknown preset {preset!r}")
        if preset and PRESETS[preset] != kind:
            raise ValueError(f"{location}: preset {preset} is incompatible with type={kind}")
        state = row.get("bam_state") or "raw"
        if state not in {"raw", "selected"}:
            raise ValueError(f"{location}: bam_state must be raw or selected")
        if kind != "duplex" and row.get("bam_state"):
            raise ValueError(f"{location}: bam_state only applies to duplex samples")
        if (row.get("fwd_pattern") and not row.get("rev_pattern")) or (row.get("rev_pattern") and not row.get("fwd_pattern")):
            raise ValueError(f"{location}: provide both orientation patterns")
        default_source, default_fwd, default_rev = "sam_flags", "forward", "reverse"
        if preset in {"UDSeq", "NanoSeq"}:
            default_fwd, default_rev = "F1R2", "F2R1"
        elif preset == "HiDEF-seq":
            default_source, default_fwd, default_rev = "read_name", r"^[^/]+/[^/]+/[^/]+/fwd(?:/|$)", r"^[^/]+/[^/]+/[^/]+/rev(?:/|$)"
        source = row.get("strand_source") or default_source
        fwd = row.get("fwd_pattern") or default_fwd
        rev = row.get("rev_pattern") or default_rev
        _validate_orientation(source, fwd, rev, location)
        family_tag = row.get("family_tag") or ""
        if family_tag and not re.fullmatch(r"[A-Za-z][A-Za-z0-9]", family_tag):
            raise ValueError(f"{location}: family_tag must be a two-character BAM tag")
        if kind != "duplex" and family_tag:
            raise ValueError(f"{location}: family_tag only applies to custom raw duplex samples")
        bam = None if row["bam"] == "." else (path.parent / row["bam"]).resolve()
        if bam is not None:
            if bam in bam_paths:
                raise ValueError(f"{location}: the same BAM is assigned to more than one sample")
            bam_paths.add(bam)
            if kind == "duplex" and state == "raw":
                if preset:
                    if family_tag:
                        raise ValueError(f"{location}: family_tag would be ignored by preset {preset}")
                    if (source, fwd, rev) != (default_source, default_fwd, default_rev):
                        raise ValueError(f"{location}: raw {preset} requires its preset orientation definitions")
                elif not family_tag or not row.get("strand_source") or not row.get("fwd_pattern"):
                    raise ValueError(f"{location}: custom raw duplex needs family_tag and explicit strand definitions")
            elif family_tag:
                raise ValueError(f"{location}: family_tag would be ignored for selected input")
        samples[sid] = Sample(sid, method, kind, bam, preset, state, source, fwd, rev, family_tag)
    if not samples or not any(sample.bam for sample in samples.values()):
        raise ValueError("samples.tsv must contain at least one query BAM")
    return samples


def load_variants(path, samples):
    required = {"chrom", "pos", "ref", "alt", "source_sample_id"}
    origins, seen = [], set()
    for number, row in _rows(path, required, required | {"callset"}):
        location = f"{path}:{number}"
        if not re.fullmatch(r"[1-9][0-9]*", row["pos"]):
            raise ValueError(f"{location}: pos must be a positive 1-based integer")
        ref, alt = row["ref"], row["alt"]
        if ref not in ("A", "C", "G", "T") or alt not in ("A", "C", "G", "T") or ref == alt:
            raise ValueError(f"{location}: requires distinct single A/C/G/T REF and ALT bases")
        if row["source_sample_id"] not in samples:
            raise ValueError(f"{location}: unknown source_sample_id {row['source_sample_id']!r}")
        allele = Allele(row["chrom"], int(row["pos"]), ref, alt)
        origin = Origin(allele, row["source_sample_id"], row.get("callset", ""))
        if origin not in seen:
            seen.add(origin)
            origins.append(origin)
    if not origins:
        raise ValueError("variants.tsv contains no candidate SNVs")
    return origins


def validate_alignments(reference, samples, alleles):
    reference = Path(reference)
    if not reference.is_file() or not Path(str(reference) + ".fai").is_file():
        raise ValueError(f"Reference FASTA and neighboring .fai are required: {reference}")
    with pysam.FastaFile(str(reference)) as fasta:
        lengths = dict(zip(fasta.references, fasta.lengths))
        for allele in alleles:
            if allele.chrom not in lengths or allele.pos > lengths[allele.chrom]:
                raise ValueError(f"Candidate outside reference: {allele}")
            actual = fasta.fetch(allele.chrom, allele.pos - 1, allele.pos).upper()
            if actual != allele.ref:
                raise ValueError(f"Reference mismatch at {allele.chrom}:{allele.pos}: table={allele.ref}, FASTA={actual}")
    chroms = {allele.chrom for allele in alleles}
    for sample in samples.values():
        if sample.bam is None:
            continue
        if not sample.bam.is_file():
            raise ValueError(f"Missing BAM for {sample.sample_id}: {sample.bam}")
        with pysam.AlignmentFile(str(sample.bam), "rb") as bam:
            if not bam.has_index():
                raise ValueError(f"Missing BAM index for {sample.sample_id}; create an index before running RIVER")
            if bam.header.to_dict().get("HD", {}).get("SO") != "coordinate":
                raise ValueError(f"BAM must declare coordinate sort order: {sample.sample_id}")
            bam_lengths = dict(zip(bam.references, bam.lengths))
            for chrom in chroms:
                if chrom not in bam_lengths:
                    raise ValueError(f"Candidate contig {chrom} missing from BAM {sample.sample_id}")
                if bam_lengths[chrom] != lengths[chrom]:
                    raise ValueError(f"Reference/BAM contig length mismatch: {sample.sample_id}, {chrom}")
    return lengths
