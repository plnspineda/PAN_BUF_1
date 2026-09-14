#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=72:00:00
#SBATCH --mem=64GB
#SBATCH --array=1-23
#SBATCH --ntasks-per-core=1

date
set -euo pipefail

module purge
module use /apps/modules/all
module load Singularity/3.10.5

samp="$1"
echo "sample_name: $samp"
mkdir -p "$samp"

# Singularity image
sif="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v3/tools/vg_v1.64.0.sif"

# Chromosome ID from array
chr=${SLURM_ARRAY_TASK_ID}
gbz_path=$(ls /hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v3/vg_index/vg_index_lloyd/chr_${chr}.fa.*.smooth.final.gbz | head -n 1)
prefix="${gbz_path%.gbz}"

echo "chr: $chr"
echo "gbz path: $gbz_path"
echo "index prefix: $prefix"

# Input reads
r1=/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v3/buffalo_shortreads/downsample_10x/"$samp"_10x_val_1.fq.gz
r2=/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v3/buffalo_shortreads/downsample_10x/"$samp"_10x_val_2.fq.gz

echo "read1: $r1"
echo "read2: $r2"

# Output file names (per chromosome)
gam="$samp/${samp}_${chr}.gam"
stats="$samp/${samp}_${chr}.stats"
pack="$samp/${samp}_${chr}"

echo "output files: $gam, $stats, $pack"

# Index files
gbz="$prefix".gbz
dist="$prefix".dist
min="$prefix".min
vg="$prefix".mod.vg

echo "gbz = $gbz"
echo "dist = $dist"
echo "min = $min"
echo "vg = $vg"

# Temporary directory
export TMPDIR=$PWD
echo "temporary directory: $TMPDIR"

# STEP: Mapping
echo "STEP: VG GIRAFFE (chr${chr})"
singularity exec "$sif" vg giraffe \
  -p \
  -t 8 \
  -Z "$gbz" \
  -d "$dist" \
  -m "$min" \
  -f "$r1" -f "$r2" > "$gam"

# STEP: Stats
echo "STEP: VG STATS chr${chr}"
singularity exec "$sif" vg stats -a "$gam" > "$stats"

# STEP: Pack
echo "STEP: VG PACK chr${chr} min mapq 10"
singularity exec "$sif" vg pack \
  -x "$vg" \
  -g "$gam" \
  --min-mapq 10 \
  -t 8 \
  -o "$pack".q10.pack

echo "STEP: VG PACK chr${chr} no filter"
singularity exec "$sif" vg pack \
  -x "$vg" \
  -g "$gam" \
  -t 8 \
  -o "$pack".pack

declare -A chr_map=(
  [1]="9,4"
  [2]="1"
  [3]="2"
  [4]="3"
  [5]="5"
  [6]="6"
  [7]="7"
  [8]="8"
  [9]="10"
  [10]="11"
  [11]="12"
  [12]="13"
  [13]="14"
  [14]="15"
  [15]="16"
  [16]="17"
  [17]="18"
  [18]="19"
  [19]="20"
  [20]="21"
  [21]="22"
  [22]="23"
  [23]="24"
)

mapped=${chr_map[$chr]}

p_args=()
IFS=',' read -r -a mapped_arr <<< "$mapped"
for m in "${mapped_arr[@]}"; do
  p_args+=(-p "ARVUO#0#${m}#0")
done

echo "Mapped -p arguments for chr${chr}:"
echo "${p_args[@]}"

# ----------------------------
# VG CALL
# ----------------------------
echo "STEP: VG CALL"
echo "singularity exec $sif vg call $vg -k $pack.pack --sample $samp ${p_args[*]} --genotype-snarls --all-snarls -t 8 > $samp/${samp}_${chr}.all.vcf"
singularity exec "$sif" vg call "$vg" -k "$pack.pack" --sample "$samp" "${p_args[@]}" --genotype-snarls --all-snarls -t 8 > "$samp/${samp}_${chr}.all.vcf"

echo "singularity exec $sif vg call $vg -k $pack.q10.pack --sample $samp ${p_args[*]} --genotype-snarls --all-snarls -t 8 > $samp/${samp}_${chr}.q10.all.vcf"
singularity exec "$sif" vg call "$vg" -k "$pack.q10.pack" --sample "$samp" "${p_args[@]}" --genotype-snarls --all-snarls -t 8 > "$samp/${samp}_${chr}.q10.all.vcf"

if ls "$gam" >/dev/null 2>&1; then
  du -sh "$gam" >> "$samp"/size.txt
  rm "$gam"
fi

echo "Done!"
date
