#!/usr/bin/env bash
set -euo pipefail #Stop immediately on error

# Usage: bash overlap_sbs.sh <outMap.vcf.gz> <sbsMap.vcf.gz> <output.vcf.gz>

OUTMAP=$1
SBSMAP=$2
OUTPUT=$3
THREADS=${4:-10}


bcftools isec -n=2 -w1 -Oz --threads "$THREADS" \
  -o "$OUTPUT" \
  "$OUTMAP" "$SBSMAP"

tabix -f -p vcf "$OUTPUT" #indexing
