# Graph-based vs linear SNP calling

To determine the r-dosage:

    bcftools stats -S - "$vcf1" "$vcf2" > "concordance.txt"

Then find the concordance part result.