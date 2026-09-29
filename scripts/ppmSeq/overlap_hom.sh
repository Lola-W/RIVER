#!/usr/bin/env bash
set -euo pipefail

# Usage: bash overlap_hom.sh <outMap.vcf.gz> <homMap.vcf.gz> <output.vcf.gz>

OUTMAP=$1
HOMMAP=$2
OUTPUT=$3
THREADS=${4:-10}

bcftools isec -n=2 -w1 -Oz --threads "$THREADS" \
  -o "$OUTPUT" \
  "$OUTMAP" "$HOMMAP"

tabix -f -p vcf "$OUTPUT"
