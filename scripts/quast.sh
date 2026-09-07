#!/bin/bash

set -euo pipefail

module purge
module use /apps/modules/all
module load quast/4.5-foss-2016uofa-Python-2.7.11

for i in genomes/*_ungapped.fa; do
  echo "$i"
  base=$(basename "$i" .fa)
  quast.py "$i" -o quast_"${base}"
done