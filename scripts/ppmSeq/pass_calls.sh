#!/usr/bin/env bash
set -euo pipefail

# Usage: bash pass_calls.sh <outMap.vcf.gz> <output.vcf.gz>

IN=$1 #input: outMap.vcf.gz
OUT=$2 #output: PASS-only variants

bcftools view -f PASS "$IN" -Oz -o "$OUT"

tabix -p vcf "$OUT"
