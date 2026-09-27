#!/usr/bin/env bash
#SBATCH --job-name=river
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH --output=river-%j.log

set -euo pipefail
if [[ $# -ne 4 ]]; then
    echo "Usage: sbatch scripts/run_river.sbatch variants.tsv samples.tsv reference.fa NEW_OUTPUT_DIR" >&2
    exit 2
fi
cd "${SLURM_SUBMIT_DIR:?Submit with sbatch from the repository root}"
export PYTHONPATH="${PWD}/src${PYTHONPATH:+:${PYTHONPATH}}"
"${RIVER_PYTHON:-python}" -m river run \
    --variants "$1" --samples "$2" --reference "$3" --outdir "$4" \
    --threads "${SLURM_CPUS_PER_TASK:-1}" --samtools "${RIVER_SAMTOOLS:-samtools}"
