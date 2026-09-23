# Graph-based vs linear SNP calling

To determine the r-dosage and weighted genotype:

    bcftools stats -S - "$vcf1" "$vcf2" > "concordance.txt"