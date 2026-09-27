"""Checked SAMtools SNV pileups and binomial Wilson intervals."""
import math
import subprocess
import tempfile
from collections import defaultdict

from .selection import observation_orientation

MPILEUP_OPTIONS = ["-Q", "13", "-q", "0", "-A", "-B", "--ff", "SECONDARY,QCFAIL,DUP,SUPPLEMENTARY", "-d", "1000000"]
WILSON_Z = 1.959963984540054
SUPPORT_THRESHOLD = 0.001


def parse_read_bases(text, ref):
    """Decode one mpileup bases field into (base-or-None, reverse) observations.

    Indel annotations are skipped without discarding their preceding SNV base.
    Deletion placeholders and reference skips retain observation slots, allowing
    associated FLAG/QNAME/tag vectors to remain aligned.
    """
    observations = []
    cursor = 0
    can_have_indel = False
    while cursor < len(text):
        char = text[cursor]
        if char == "^":
            if cursor + 2 >= len(text):
                raise ValueError("Malformed mpileup read-start marker")
            cursor += 2
            can_have_indel = False
            continue
        if char == "$":
            if not can_have_indel:
                raise ValueError("Malformed mpileup read-end marker")
            cursor += 1
            can_have_indel = False
            continue
        if char in "+-":
            if not can_have_indel:
                raise ValueError("Indel annotation has no preceding observation")
            start = cursor + 1
            end = start
            while end < len(text) and text[end].isdigit():
                end += 1
            if end == start:
                raise ValueError("Indel annotation lacks a length")
            length = int(text[start:end])
            if length <= 0 or end + length > len(text):
                raise ValueError("Truncated or empty mpileup indel annotation")
            payload = text[end:end + length]
            if any(base not in "ACGTNacgtn*#" for base in payload):
                raise ValueError("Invalid base in mpileup indel annotation")
            cursor = end + length
            continue
        if char == ".":
            observations.append((ref.upper(), False))
        elif char == ",":
            observations.append((ref.upper(), True))
        elif char in "ACGTNacgtn":
            observations.append((char.upper(), char.islower()))
        elif char in "*#<>":
            observations.append((None, char in "#<"))
        else:
            raise ValueError(f"Unsupported mpileup character: {char!r}")
        cursor += 1
        can_have_indel = True
    return observations


def wilson_interval(alt_count, depth):
    if not isinstance(alt_count, int) or not isinstance(depth, int) or depth < 0 or not 0 <= alt_count <= depth:
        raise ValueError("Wilson counts must be integers with 0 <= ALT <= depth")
    if depth == 0:
        return None, None, None
    proportion = alt_count / depth
    squared = WILSON_Z * WILSON_Z
    denominator = 1.0 + squared / depth
    center = proportion + squared / (2.0 * depth)
    offset = WILSON_Z * math.sqrt(proportion * (1.0 - proportion) / depth + squared / (4.0 * depth * depth))
    return proportion, max(0.0, (center - offset) / denominator), min(1.0, (center + offset) / denominator)


def _vectors(field, count, description):
    values = field.split(",")
    if len(values) != count:
        raise ValueError(f"mpileup {description} vector length {len(values)} differs from observation count {count}; commas in read names/tag values are unsupported")
    return values


def parse_pileup_line(line, sample):
    """Return (chrom, pos, reference, counts[base][orientation])."""
    fields = line.rstrip("\n").split("\t")
    tag = sample.strand_source.startswith("tag:")
    if len(fields) != (9 if tag else 8):
        raise ValueError(f"Expected {9 if tag else 8} mpileup columns, found {len(fields)}")
    chrom, position, ref, depth = fields[:4]
    position, depth = int(position), int(depth)
    if depth < 0 or position <= 0 or ref.upper() not in "ACGTN" or len(ref) != 1:
        raise ValueError("Invalid mpileup locus/reference/depth")
    if depth == 0:
        if fields[4] != "*" or fields[5] != "*":
            raise ValueError("Unexpected observations in zero-depth mpileup row")
        return chrom, position, ref.upper(), {}
    observations = parse_read_bases(fields[4], ref)
    if len(observations) != depth or len(fields[5]) != depth:
        raise ValueError(f"mpileup depth/observation/quality lengths differ at {chrom}:{position}")
    names = _vectors(fields[6], depth, "QNAME")
    flags = _vectors(fields[7], depth, "FLAG")
    values = _vectors(fields[8], depth, "strand tag") if tag else [None] * depth
    counts = defaultdict(lambda: [0, 0, 0])
    for (base, _), name, flag, value in zip(observations, names, flags, values):
        orientation = observation_orientation(sample, int(flag), name, None if value == "*" else value)
        if base is not None:
            counts[base][orientation] += 1
    return chrom, position, ref.upper(), dict(counts)


def _row(allele, sample, counts, status):
    ref = counts.get(allele.ref, [0, 0, 0])
    alt = counts.get(allele.alt, [0, 0, 0])
    nr, na = sum(ref), sum(alt)
    vaf, lower, upper = wilson_interval(na, nr + na)
    supported = na >= 1 if sample.type == "single_cell" else lower is not None and lower > SUPPORT_THRESHOLD
    return {"chrom": allele.chrom, "pos": allele.pos, "ref": allele.ref, "alt": allele.alt,
            "sample_id": sample.sample_id, "method": sample.method, "type": sample.type,
            "ref_count": nr, "alt_count": na, "depth": nr + na, "vaf": vaf,
            "lower_ci": lower, "upper_ci": upper, "qualifying_support": supported,
            "query_status": status, "coverage_status": "covered" if nr + na else "no_ref_alt_depth",
            "ref_forward": ref[1], "ref_reverse": ref[2], "ref_unclassified": ref[0],
            "alt_forward": alt[1], "alt_reverse": alt[2], "alt_unclassified": alt[0]}


def query_sample(sample, bam, reference, alleles, bed_path, samtools, commands):
    """Stream an indexed chromosome query once per chromosome with candidate BED."""
    by_chrom = defaultdict(dict)
    for allele in alleles:
        by_chrom[allele.chrom].setdefault(allele.pos, []).append(allele)
    for chrom, positions in by_chrom.items():
        extra = "QNAME,FLAG"
        if sample.strand_source.startswith("tag:"):
            extra += "," + sample.strand_source[4:]
        region = f"{chrom}:{min(positions)}-{max(positions)}"
        command = [samtools, "mpileup", "-r", region, "-l", str(bed_path), "-f", str(reference)] + MPILEUP_OPTIONS + ["--output-extra", extra, str(bam)]
        commands.append(command)
        seen = set()
        with tempfile.TemporaryFile(mode="w+") as errors:
            process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=errors, text=True)
            try:
                for line in process.stdout:
                    found_chrom, pos, reference_base, counts = parse_pileup_line(line, sample)
                    if found_chrom != chrom or pos not in positions:
                        raise ValueError(f"Unexpected mpileup locus {found_chrom}:{pos}")
                    if pos in seen:
                        raise ValueError(f"Duplicate mpileup locus {chrom}:{pos}")
                    seen.add(pos)
                    for allele in positions[pos]:
                        if reference_base != allele.ref:
                            raise ValueError(f"mpileup/reference mismatch at {chrom}:{pos}")
                        yield _row(allele, sample, counts, "queried")
                returncode = process.wait()
                if returncode:
                    errors.seek(0)
                    raise RuntimeError(f"SAMtools mpileup failed for {sample.sample_id}, {chrom} (exit {returncode}): {errors.read().strip()}")
            finally:
                if process.poll() is None:
                    process.terminate()
                    process.wait()
                process.stdout.close()
        for pos, site_alleles in positions.items():
            if pos not in seen:
                for allele in site_alleles:
                    yield _row(allele, sample, {}, "no_pileup_observations")
