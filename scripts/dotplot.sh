#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 32
#SBATCH --time=24:00:00
#SBATCH --mem=64GB
#SBATCH --array=0-24
#SBATCH --ntasks-per-core=1

date
module purge
module use /apps/modules/all
module load Java/17.0.6

refdir="buffalo_pangenome/input-files/assemblies_per_chr"
qrys=(buffalo_pangenome/input-files/assemblies_per_chr/RVNR1*.fna)
qry=${qrys[$SLURM_ARRAY_TASK_ID]}
base=$(basename "$qry" .fna)
chr=$(echo "$base" | sed 's/RVNR1//')

echo "$qry"

java -Xmx16g -cp tools/gepard/dist/Gepard-2.1.jar org.gepard.client.cmdline.CommandLine \
  -seq "$refdir"/RVAZ1"$chr".fna "$qry" \
  -matrix tools/gepard/resources/matrices/edna.mat \
  -outfile "$base".png -word 100 -window 10 -lower 10

date