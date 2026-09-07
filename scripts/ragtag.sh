#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 24
#SBATCH --time=24:00:00
#SBATCH --mem=64GB
#SBATCH --array=0,2-3
#SBATCH --ntasks-per-core=1

date
#conda activate ragtag
module purge
module use /apps/modules/all
module load networkx/2.6.3-foss-2021b
module load minimap2/2.26-GCCcore-11.2.0
module load MUMmer/4.0.0beta-foss-2016b
module load numpy/1.18.1-foss-2016b-Python-3.7.0
module load Pysam/0.15.4-foss-2016b-Python-3.7.0
module load intervaltree/0.1-GCCcore-11.2.0
module load SAMtools/1.17-GCC-11.2.0

export PATH=/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/tools/RagTag/:$PATH

ref="REFERENCES/UOA_WB_1/UOA_WB_1.chr.fa"
files=(buffalo_pangenome/input-files/assemblies_genome/fazeela/assembly_out/*hap*)
qry=${files[$SLURM_ARRAY_TASK_ID]}
echo "Processing $qry..."

base=$(basename "$qry" .p_ctg.fa)

mkdir -p "$base" && cd "$base" || exit
ragtag.py scaffold "$ref" "$qry"

echo "Done."
date