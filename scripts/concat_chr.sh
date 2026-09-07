#!/bin/bash

input_dir="$1"

if [[ -z "$input_dir" || ! -d "$input_dir" ]]; then
  echo "Provide input directory containing the chromosomes"
  echo "Usage: $0 <input_dir>"
  exit 1
fi

echo "$input_dir"

sw=({2..23})
rv=(1 2 3 5 6 7 8 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24)
rv_asm=( RVND1 RVND2 RVAZ1 RVAZ2 RVNR1 RVNR2 ARVUO0 )
sw_asm=( SWPC1 SWPC2 )

for ((idx=0; idx<${#sw[@]}; idx++)); do
    rv_out=""
    sw_out=""

    for ((i=0; i<${#rv_asm[@]}; i++)); do
        if [[ "${rv_asm[i]}" == "ARVUO0" ]]; then
            rv_out+=" ${input_dir}/${rv_asm[i]}${rv[idx]}.fa"
        else
            rv_out+=" ${input_dir}/${rv_asm[i]}${rv[idx]}_ungapped.fa"
        fi
    done

    for ((i=0; i<${#sw_asm[@]}; i++)); do
        sw_out+=" ${input_dir}/${sw_asm[i]}${sw[idx]}_ungapped.fa"
    done

    seqs=$(echo "$rv_out $sw_out")
    echo "$seqs > chr_${sw[idx]}.fa"
    cat $seqs > chr_"${sw[idx]}".fa
done

## for chromosome X
dir=buffalo_pangenome/input-files/assemblies_per_chr3/X_Y
dir_sw=buffalo_pangenome_v2/input_files/haplotypes

cat "$dir"/RVND2X_ungapped.fa "$dir"/RVND1X_ungapped.fa \
ARVUO0X.fa \
"$dir"/RVAZ1X_ungapped.fa "$dir"/RVAZ2X_ungapped.fa \
"$dir"/RVNR1X_ungapped.fa "$dir"/RVNR2X_ungapped.fa \
"$dir_sw"/X_1_ungapped.fa "$dir_sw"/X_2_ungapped.fa  > chr_X.fa

## for chromosome 1
cat "$dir"/RVND24_ungapped.fa "$dir"/RVND29_ungapped.fa "$dir"/RVND14_ungapped.fa "$dir"/RVND19_ungapped.fa ARVUO04.fa ARVUO09.fa "$dir"/RVAZ14_ungapped.fa "$dir"/RVAZ19_ungapped.fa "$dir"/RVAZ24_ungapped.fa "$dir"/RVAZ29_ungapped.fa "$dir"/RVNR14_ungapped.fa "$dir"/RVNR19_ungapped.fa "$dir"/RVNR24_ungapped.fa "$dir"/RVNR29_ungapped.fa SWPC11_ungapped.fa SWPC21_ungapped.fa > chr_1.fa

