#!/bin/bash
#SBATCH -p icelake
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --time=72:00:00
#SBATCH --mem=64GB

set -euo pipefail
date

module purge
module use /apps/modules/all
module load GATK/4.6.1.0-GCCcore-12.3.0-Java-17.0.6

INPUT_VCF="combined_snps.vcf.gz"
OUTPUT_VCF="combined_snps.filter.vcf.gz"

echo "input: $INPUT_VCF"
echo "output: $OUTPUT_VCF"

gatk VariantFiltration \
    -V ${INPUT_VCF} \
    --filter-expression "QD < 2.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
    --filter-name "filter1" \
    --filter-expression "DP > 400" \
    --filter-name "highDP" \
    --filter-expression "DP < 140" \
    --filter-name "lowDP" \
    --genotype-filter-expression "GQ < 20.0" \
    --genotype-filter-name "lowGQ" \
    --cluster-size 3 \
    --cluster-window-size 10 \
    -O ${OUTPUT_VCF}

gatk --java-options "-Xmx48g" IndexFeatureFile -I "${OUTPUT_VCF}"

echo "Done!"
date

