#!/usr/bin/env bash
set -euo pipefail
if [[ $# -gt 1 ]]; then
    echo "Usage: bash examples/run_demo.sh [NEW_DEMO_DIRECTORY]" >&2
    exit 2
fi
example_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
demo_dir="${1:-${example_dir}/demo}"
if [[ -e "$demo_dir" ]]; then
    echo "Demo directory already exists; choose a new directory: $demo_dir" >&2
    exit 1
fi
python="${RIVER_PYTHON:-python}"
"$python" "$example_dir/generate_synthetic.py" --outdir "$demo_dir/data"
"$python" -m river run \
    --variants "$demo_dir/data/variants.tsv" \
    --samples "$demo_dir/data/samples.tsv" \
    --reference "$demo_dir/data/reference.fa" \
    --outdir "$demo_dir/results" \
    --threads 1 --samtools "${RIVER_SAMTOOLS:-samtools}"
"$python" "$example_dir/check_demo.py" "$demo_dir/results"
