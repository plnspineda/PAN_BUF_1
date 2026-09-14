#!/bin/bash -l
#SBATCH -p a100cpu
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --time=01:00:00
#SBATCH --mem=12GB
#SBATCH --array=0-3  # Replace XX with number of lines - 1 in fasta_list.txt

# Read the line corresponding to the array task ID
#LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" list_fasta_for_SVIM.txt)
LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" list2_fasta_for_SVIM.txt)
#LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID + 1))p" list3_fasta_for_SVIM.txt) #did not use because I used copy_vcfs.sh

# Notes: remember to samtools index *.fa

# Parse columns
ID=$(echo "$LINE" | cut -f1)
hap1=$(echo "$LINE" | cut -f2)
hap2=$(echo "$LINE" | cut -f3)

echo $ID
echo $hap1
echo $hap2

# Run the pipeline script
# --time=48:00:00 --mem=96GB
#./SVIM_asm.sh "$ID" "$hap1" "$hap2"

# below is just for copying the vcf and put them in one place
# --time=01:00:00 --mem=12GB
./SVIM_asm_copy_vcf.sh "$ID" "$hap1" "$hap2"
