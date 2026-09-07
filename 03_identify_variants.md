# Identifying variants of the pangenome

Github page: https://github.com/vgteam/vg

`vg toolkit` is a multifunction tool for working with pangenome graphs.

1. Concatenate the vcf files using `scripts/bcf_concat.sh`.
2. Filter and clean the vcf file using vcflib tools `scripts/vcfbub.sh`
3. Separate the variants using `scripts/seperate_variants_v4.sh`

Get vcf stats for any vcf files using `scripts/bcf_get_stats.sh`

Can also classify the different variants using `scripts/classify_bubbles_wLengthAlleleFreq.sh` and `scripts/allele_metrics.sh`

These are the different variant classifications:

| Allelic variation | Variation type | Classification | Definition |
|---|---|---|---|
| biallelic | SNP | SNP | Single nucleotide polymorphism with only 1-bp allele for both REF and ALT. |
| biallelic | SNP | SNP other | Single nucleotide polymorphism with longer than 1-bp allele for either REF and ALT. For example, an allele can be (REF) CAAT and (ALT) CCAT, which is a SNP but the alleles are longer than 1-bp. |
| biallelic | MNP | MNP | Multinucleotide polymorphism with 2-bp allele for both REF and ALT. |
| biallelic | InDel | small InDel | Insertion or deletion variants with alleles shorter than 50-bp for both REF and ALT. |
| biallelic | SCV | SCV | Small complex variants that are not InDel and shorter than 50-bp for both the REF and ALT allele. |
| biallelic | SV | SV InDel | Structural variants that are InDel and longer than 50-bp for either the REF and ALT allele. |
| biallelic | SV | SV complex | Structural variants that are not InDel and longer than 50-bp for either the REF and ALT allele. |
| multiallelic | SNP | SNP | Single nucleotide polymorphism with only 1-bp allele for both REF and ALT, and contains more than one ALT alleles. |
| multiallelic | SNP | SNP other | Single nucleotide polymorphism with longer than 1-bp allele for either REF and ALT, and contains more than one ALT alleles. For example, an allele can be (REF) CAAT and (ALT) CCAT, which is a SNP but the alleles are longer than 1-bp. |
| multiallelic | MNP | MNP | Multinucleotide polymorphism with 2-bp allele for both REF and ALT, and contains more than one ALT alleles. |
| multiallelic | InDel | small InDel | Insertion or deletion variants with alleles shorter than 50-bp for both REF and ALT, and contains more than one ALT alleles. |
| multiallelic | SCV | SCV | Small complex variants that are not InDel and shorter than 50-bp for both the REF and ALT allele, and contains more than one ALT alleles. |
| multiallelic | SV | SV InDel | Structural variants that are InDel and longer than 50-bp for either the REF and ALT allele, and contains more than one ALT alleles. |
| multiallelic | SV | SV complex | Structural variants that are not InDel and longer than 50-bp for either the REF and ALT allele, and contains more than one ALT alleles. |