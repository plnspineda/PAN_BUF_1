#!/bin/bash

Ref=/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v3/reference_fa/ARVUO_withX.fa
OUT=/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v3/SVIM/OUT

AssID=$1
Ass1=$2
Ass2=$3

mkdir -p $OUT/$AssID/
mkdir -p $OUT/$AssID/SVIM

module use /apps/modules/all
module load Anaconda3/2023.03
module load SAMtools/1.17-GCC-11.2.0

eval "$(conda shell.bash hook)"

conda activate svimasm_env

minimap2 -a -x asm5 --cs -r2k -t 16 $Ref $Ass1 | \
samtools sort -@ 16 -o "$OUT/$AssID/SVIM/${AssID}.hap1.sorted.bam"

minimap2 -a -x asm5 --cs -r2k -t 16 $Ref $Ass2 | \
samtools sort -@ 16 -o "$OUT/$AssID/SVIM/${AssID}.hap2.sorted.bam"

samtools index "$OUT/$AssID/SVIM/${AssID}.hap1.sorted.bam"
samtools index "$OUT/$AssID/SVIM/${AssID}.hap2.sorted.bam"

svim-asm diploid --sample "$AssID" --min_sv_size 50 --query_names $OUT/$AssID/SVIM \
  "$OUT/$AssID/SVIM/${AssID}.hap1.sorted.bam" \
  "$OUT/$AssID/SVIM/${AssID}.hap2.sorted.bam" \
  $Ref

conda deactivate

echo $AssID Finished!
