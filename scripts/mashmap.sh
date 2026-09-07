#!/bin/bash

module purge
module use /apps/modules/all
module load mashmap/2.0
module load gnuplot/5.4.2-GCCcore-11.2.0

outdir="mashmap_checkaln"
ref="SWPC.fna"
files=(./*fna)
qry=${files[$SLURM_ARRAY_TASK_ID]}

echo -e "Processing ${qry} as query with $ref as reference genome."

baseref=$(basename "$ref" .fna)
baseqry=$(basename ${qry} .fna)

mkdir -p "$outdir"/"$baseqry"-"$baseref"
cd "$outdir"/"$baseqry"-"$baseref"
mashmap -r $ref -q $qry --pi 95 -t 4 -f one-to-one -s 5000000 -o "$baseref"_"$baseqry"_mashmap.out
tools/MashMap/scripts/generateDotPlot png medium "$baseref"_"$baseqry"_mashmap.out

awk '{if ($1 != $6) print "\033[93mWarning:\033[0m Chromosome "$1" and "$6" are not equal"}' "$baseref"_"$baseqry"_mashmap.out

echo -e "\n Finished ${qry}."