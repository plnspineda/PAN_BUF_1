#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 12
#SBATCH --time=12:00:00
#SBATCH --mem=128GB

set -euo pipefail

date
module purge
module use /apps/modules/all
module load Singularity/3.10.5
module load SAMtools/1.17-GCC-11.2.0

chr=X
gfa_dir=/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v2/pggb_pangenie/haplotyped_2025.08.14/chr_${chr}_p95_s10k
gfa=${gfa_dir}/chr_${chr}*.gfa
echo "$gfa"

for f in $gfa; do
    base_gfa=$(basename "$f" .gfa)
    echo "$base_gfa"
done

### index graph
echo "→ vg convert"
singularity exec /hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/vg_v1.64.0/vg_v1.64.0.sif vg \
 convert -g $gfa > "$base_gfa".vg

echo "→ vg mod"
singularity exec /hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/vg_v1.64.0/vg_v1.64.0.sif vg \
 mod -X 256 "$base_gfa".vg > "$base_gfa".mod.vg

echo "→ vg prune"
singularity exec /hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/vg_v1.64.0/vg_v1.64.0.sif vg \
 prune "$base_gfa".mod.vg -t 8 > "$base_gfa".pruned.vg

echo "→ vg index gbwt"
mkdir -p tmp_"$base_gfa"
singularity exec /hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/vg_v1.64.0/vg_v1.64.0.sif vg \
 index -x "$base_gfa".xg -g "$base_gfa".gbwt "$base_gfa".pruned.vg -t 8 -b tmp_"$base_gfa"

echo "→ vg index dist"
singularity exec /hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/vg_v1.64.0/vg_v1.64.0.sif vg \
 index -j "$base_gfa".dist -t 8 "$base_gfa".pruned.vg

echo "→ vg gbwt"
singularity exec /hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/vg_v1.64.0/vg_v1.64.0.sif vg \
 gbwt -x "$base_gfa".xg -o "$base_gfa".gbwt --path-cover --num-threads 8

echo "→ vg gbwt gbz"
singularity exec /hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/vg_v1.64.0/vg_v1.64.0.sif vg \
 gbwt -x "$base_gfa".xg -g  "$base_gfa".gbz --gbz-format "$base_gfa".gbwt

echo "→ vg minimizer"
singularity exec /hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome/tools/vg_v1.64.0/vg_v1.64.0.sif vg \
 minimizer -d "$base_gfa".dist --gbwt-name "$base_gfa".gbwt -o "$base_gfa".min "$base_gfa".mod.vg -t 8

rm -r  tmp_"$base_gfa"

echo "Done!"

