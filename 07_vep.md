# Prediction of high impact variants


## orthofinder

    conda activate orthofinder
    orthofinder -f genomes/ -t 24

## variant effect predictor

    #!/bin/bash
    date
    module purge
    module use /apps/modules/all
    module load Singularity/3.10.5

    vep="/hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/buffalo_pangenome_v3/tools/ensembl-vep.sif"
    seq="/hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/REFERENCES/UOA_WB_1/UOA_WB_1.fa.gz"
    gene="/hpcfs/groups/phoenix-hpc-avsci/Davies_Informatics/REFERENCES/UOA_WB_1/UOA_WB_1.sorted.gff.gz"
    vcf="$1" ## contains the variants for vep

    if [[ -z "$vcf" ]]; then
        echo "Error: No VCF file provided."
        echo "Usage: sbatch $0 <input.vcf.gz>"
        exit 1
    fi

    echo "vcf file: $vcf"
    echo "sequence: $seq"
    echo "gene: $gene"

    base="$(basename "$vcf" .vcf.gz)"
    dir=vep_ncbi_"$base"

    gunzip -c "$vcf" | sed 's/ARVUO#0#//g' > tmp_"$base"

    echo "Running..."
    mkdir -p "$dir"
    singularity exec $vep vep --fork 12 -i tmp_"$base" -o "$dir"/"$base".res.txt --gff $gene --fasta $seq

    rm tmp_"$base"
    echo "Done!"
    date

Combine output data for all variant types using `scripts/combine_all_data_whole.R`

Combine all impact variants using `scripts/combine_all_impact_2.sh`

Plot the graphs using `plot_scripts/VEP_analyse_impactR.R` and `plot_scripts/plot_bargraph_VEP.R`