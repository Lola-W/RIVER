#!/usr/bin/env python3
"""Compare all demo output columns with the checked-in expected tables."""
import argparse
import csv
import json
import math
from pathlib import Path


FLOAT_FIELDS = {"vaf", "lower_ci", "upper_ci", "proportion_support_0",
                "proportion_support_1", "proportion_support_ge2"}


def read_table(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames, list(reader)


def check(results):
    expected = Path(__file__).with_name("expected")
    for name in ("read_support.tsv", "variant_support.tsv", "summary.tsv"):
        fields, wanted = read_table(expected / name)
        actual_fields, actual = read_table(results / name)
        if actual_fields != fields or len(actual) != len(wanted):
            raise ValueError(f"{name}: header or row count differs from expected")
        for number, (observed, reference) in enumerate(zip(actual, wanted), 2):
            if None in observed or any(value is None for value in observed.values()):
                raise ValueError(f"{name}:{number}: malformed TSV row")
            for field in fields:
                got, want = observed[field], reference[field]
                matches = got == want
                if field in FLOAT_FIELDS and got != "NA" and want != "NA":
                    matches = math.isclose(float(got), float(want), rel_tol=1e-12, abs_tol=1e-12)
                if not matches:
                    raise ValueError(f"{name}:{number}:{field}: expected {want!r}, got {got!r}")
    metadata = json.loads((results / "run.json").read_text())
    if metadata.get("status") != "complete":
        raise ValueError("run.json: expected status=complete")
    prep = metadata["preparation"]["duplex_library"]
    if (prep["eligible_reads"], prep["written_reads"]) != (21, 20):
        raise ValueError("run.json: expected ppmSeq selection to retain 20 of 21 reads")
    print("Demo matches expected tables: 12 read rows, 4 variant rows, 2 summary rows.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results", type=Path)
    args = parser.parse_args()
    try:
        check(args.results)
    except (ValueError, OSError, KeyError) as exc:
        parser.exit(1, f"Demo check failed: {exc}\n")


if __name__ == "__main__":
    main()
