"""Preset and user-defined duplex read selection."""
import re
import subprocess

import pysam

EXCLUDE_FLAGS = 3840
PRIMARY_CHROMS = {str(number) for number in range(1, 23)} | {"X"}


def passes_read_qc(rec) -> bool:
    if rec.is_secondary:
        return False
    if rec.is_supplementary:
        return False

    if rec.is_qcfail:
        return False

    if not rec.is_proper_pair:
        return False

    if rec.has_tag("DT"):
        return False

    ct = rec.cigartuples
    if ct and len(ct) > 0:
        op, length = ct[0]
        if length == 4:
            return False

    return True


def canonical_umi_pair(rec) -> str:
    """Normalize the terminal UMI pair using read number and orientation."""
    bc = rec.query_name.rsplit("_", 1)[-1]
    if "+" not in bc:
        return bc
    bc1, bc2 = bc.split("+", 1)

    if (rec.is_read1 and rec.is_forward) or (rec.is_read2 and rec.is_reverse):
        return f"{bc1}+{bc2}"
    else:
        return f"{bc2}+{bc1}"


def classify_f1r2_f2r1(rec) -> int:
    """Return F1R2=1 or F2R1=2 from read number and alignment orientation."""
    if (rec.is_forward and rec.is_read1) or (rec.is_reverse and rec.is_read2):
        return 1
    return 2


def observation_orientation(sample, flag, name="", tag_value=None):
    """Return 1/2 for the configured labels, or 0 for an unclassified read."""
    if sample.strand_source == "sam_flags":
        reverse = bool(flag & 16)
        if sample.fwd_pattern == "forward":
            return 2 if reverse else 1
        read1, read2 = bool(flag & 64), bool(flag & 128)
        if read1 == read2:
            return 0
        return 1 if (read1 and not reverse) or (read2 and reverse) else 2
    value = name if sample.strand_source == "read_name" else tag_value
    if value is None:
        return 0
    forward = re.search(sample.fwd_pattern, str(value)) is not None
    reverse = re.search(sample.rev_pattern, str(value)) is not None
    if forward and reverse:
        raise ValueError(f"{sample.sample_id}: both strand patterns match {value!r}")
    return 1 if forward else 2 if reverse else 0


def _record_orientation(sample, record):
    value = None
    if sample.strand_source.startswith("tag:"):
        tag = sample.strand_source[4:]
        if record.has_tag(tag):
            value = record.get_tag(tag)
    return observation_orientation(sample, record.flag, record.query_name, value)


def _primary_contigs(bam):
    return [chrom for chrom in bam.references if chrom.removeprefix("chr") in PRIMARY_CHROMS]


def _nano_key(record, sample):
    name = record.query_name or ""
    tail = name.rsplit("_", 1)[-1]
    if "_" not in name or tail.count("+") != 1 or not all(tail.split("+")):
        raise ValueError(f"{sample.sample_id}: UDSeq/NanoSeq read name must end in _UMI1+UMI2: {name!r}")
    return record.reference_start, canonical_umi_pair(record)


def _custom_key(record, sample):
    if not record.has_tag(sample.family_tag):
        raise ValueError(f"{sample.sample_id}: eligible read lacks family tag {sample.family_tag}: {record.query_name}")
    value = record.get_tag(sample.family_tag)
    if not isinstance(value, (str, int)) or value == "":
        raise ValueError(f"{sample.sample_id}: family tag must contain a nonempty string or integer")
    return value


def prepare_bam(sample, output, samtools, threads=1):
    """Select duplex reads into a new, sorted/indexed BAM; never modify inputs."""
    output.parent.mkdir(parents=True, exist_ok=True)
    unsorted = output.with_suffix(".unsorted.bam")
    stats = {"eligible_reads": 0, "written_reads": 0, "duplex_families": 0, "unclassified_reads": 0}
    with pysam.AlignmentFile(str(sample.bam), "rb", threads=threads) as source:
        contigs = _primary_contigs(source) if sample.preset else list(source.references)
        stats["contigs"] = contigs
        if not contigs:
            raise ValueError(f"{sample.sample_id}: no chromosomes 1–22/X available for preset selection")
        with pysam.AlignmentFile(str(unsorted), "wb", header=source.header, threads=threads) as target:
            if sample.preset == "HiDEF-seq":
                masks = {}
                stats["recognized_label_reads"] = 0
                for chrom in contigs:
                    for rec in source.fetch(chrom):
                        if rec.flag & EXCLUDE_FLAGS:
                            continue
                        stats["eligible_reads"] += 1
                        fields = (rec.query_name or "").split("/")
                        if len(fields) < 2:
                            stats["unclassified_reads"] += 1
                            continue
                        label = fields[3] if len(fields) >= 4 else ""
                        bit = 1 if label == "fwd" else 2 if label == "rev" else 0
                        stats["recognized_label_reads"] += int(bit != 0)
                        stats["unclassified_reads"] += int(bit == 0)
                        masks[fields[1]] = masks.get(fields[1], 0) | bit
                if stats["eligible_reads"] and not stats["recognized_label_reads"]:
                    raise ValueError(f"{sample.sample_id}: no eligible HiDEF-seq reads have fwd/rev in slash-delimited query-name field 4; check the preset and input BAM")
                stats["duplex_families"] = sum(mask == 3 for mask in masks.values())
                for chrom in contigs:
                    for rec in source.fetch(chrom):
                        fields = (rec.query_name or "").split("/")
                        if not rec.flag & EXCLUDE_FLAGS and len(fields) >= 2 and masks.get(fields[1]) == 3:
                            target.write(rec)
                            stats["written_reads"] += 1
            elif sample.preset == "ppmSeq":
                stats["eligible_with_required_tags"] = 0
                for chrom in contigs:
                    for rec in source.fetch(chrom):
                        if rec.flag & EXCLUDE_FLAGS:
                            continue
                        stats["eligible_reads"] += 1
                        tagged = rec.has_tag("st") and rec.has_tag("et")
                        stats["eligible_with_required_tags"] += int(tagged)
                        if tagged and rec.get_tag("st") == "MIXED" and rec.get_tag("et") == "MIXED":
                            target.write(rec)
                            stats["written_reads"] += 1
                if stats["eligible_reads"] and not stats["eligible_with_required_tags"]:
                    raise ValueError(f"{sample.sample_id}: no eligible ppmSeq reads contain both st and et tags; check the preset and input BAM")
            else:
                nano = sample.preset in {"UDSeq", "NanoSeq"}
                for chrom in contigs:
                    masks = {}
                    for rec in source.fetch(chrom):
                        if (not passes_read_qc(rec)) if nano else bool(rec.flag & EXCLUDE_FLAGS):
                            continue
                        stats["eligible_reads"] += 1
                        key = _nano_key(rec, sample) if nano else _custom_key(rec, sample)
                        bit = classify_f1r2_f2r1(rec) if nano else _record_orientation(sample, rec)
                        stats["unclassified_reads"] += int(bit == 0)
                        masks[key] = masks.get(key, 0) | bit
                    stats["duplex_families"] += sum(mask == 3 for mask in masks.values())
                    for rec in source.fetch(chrom):
                        if (not passes_read_qc(rec)) if nano else bool(rec.flag & EXCLUDE_FLAGS):
                            continue
                        key = _nano_key(rec, sample) if nano else _custom_key(rec, sample)
                        if masks.get(key) == 3 and (nano or _record_orientation(sample, rec) != 0):
                            target.write(rec)
                            stats["written_reads"] += 1
                if not nano and stats["eligible_reads"] and stats["unclassified_reads"] == stats["eligible_reads"]:
                    raise ValueError(f"{sample.sample_id}: no eligible custom duplex reads match either strand definition; check strand_source and forward/reverse patterns")
    subprocess.run([samtools, "sort", "-@", str(threads), "-o", str(output), str(unsorted)], check=True)
    subprocess.run([samtools, "index", "-@", str(threads), str(output)], check=True)
    unsorted.unlink()
    return stats
