#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 16
#SBATCH --time=24:00:00
#SBATCH --mem=48GB
#SBATCH --array=1-23
#SBATCH --ntasks-per-core=1

set -euo pipefail

date
module purge
module use /apps/modules/all
module load Singularity/3.10.5
module load SAMtools/1.17-GCC-11.2.0

chr=$SLURM_ARRAY_TASK_ID
asm=/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v2/input_files/per_chr_3colnames/chr_"$chr".fa
echo "$asm"

base=$(basename "$asm" .fa)

#mkdir -p "$base" || cd "$base"
samtools faidx "$asm"
singularity run \
  -B /hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v2/tools/pggb \
  /hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v2/tools/pggb_latest.sif \
  pggb -i "$asm" -n 5 -p 95 -s 10k -t 16 -V ARVUO -o "$base"_p95_s10k

date