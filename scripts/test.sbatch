#!/usr/bin/env bash
#SBATCH --job-name=river-tests
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --output=river-tests-%j.log

set -euo pipefail
# Submit from the repository root with the appropriate account/partition.
cd "${SLURM_SUBMIT_DIR:?Submit with sbatch from the repository root}"
export PYTHONPATH="${PWD}/src${PYTHONPATH:+:${PYTHONPATH}}"
"${RIVER_PYTHON:-python}" -m unittest discover -s tests -v
