#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=24:00:00
#SBATCH --mem=64GB

date
set -euo pipefail

module purge
module use /apps/modules/all
module load BCFtools/1.17-GCC-11.2.0
module load Jellyfish/2.2.6-foss-2016b
module load GLib/2.77.1-GCCcore-12.3.0

if [[ -z "${1:-}" || -z "${2:-}" || -z "${3:-}" ]]; then
    echo "Usage: $0 <cleaned fastq directory> <sample_name> <out_dir>"
    echo "Example: sbatch $0 /home/reads sample1 swamp"
    exit 1
fi

dir="$1"
samp="$2"
out="$3"

echo "command: sbatch $0 $1 $2 $3"
echo "input fastq directory: $dir"
echo "input sample name: $samp"

vcf="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v3/pggb_concat_VCF/pggb_with_p93_chr8/ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.vcf"
ref="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v3/reference_fa/ARVUO_withX.fa"

echo "fixed input vcf: $vcf"
echo "fixed input ref: $ref"

# input reads
shopt -s nullglob

r1_gz=( "$dir"/"$samp"*val_1.fq.gz )
r2_gz=( "$dir"/"$samp"*val_2.fq.gz )

if (( ${#r1_gz[@]} != 1 || ${#r2_gz[@]} != 1 )); then
    echo "ERROR: Expected exactly one R1 and one R2 file"
    echo "R1 found: ${r1_gz[*]}"
    echo "R2 found: ${r2_gz[*]}"
    exit 1
fi

echo "R1 gz: ${r1_gz[0]}"
echo "R2 gz: ${r2_gz[0]}"

r1="${samp}_R1.fq"
r2="${samp}_R2.fq"

# decompress reads
echo "Decompressing reads..."
gunzip -c "${r1_gz[0]}" > "$r1"
gunzip -c "${r2_gz[0]}" > "$r2"

echo "Running: PanGenie"
#/usr/bin/time -v PanGenie-index -v "$vcf" -r "$ref" -t 8 -o pan
echo /usr/bin/time -v PanGenie -f pan_index/pan -i "$r1" "$r2" -j 8 -t 8 -s "$samp" -o "$out/panvar_$samp"
/usr/bin/time -v PanGenie -f pan_index/pan -i "$r1" "$r2" -j 8 -t 8 -s "$samp" -o "$out"/panvar_"$samp"

echo "Cleaning up uncompressed files..."
rm -f "$r1" "$r2"

echo "Done!"
date
