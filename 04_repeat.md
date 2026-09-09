# Repeat Analysis

The SVs identified from PGGB after vcfbub filtering were converted into FASTA format by first normalizing the VCF to separate multiallelic sites, then extracting the longest allele from either the reference or alternate alleles as the representative sequence. 

Normalise the vcf file

    bcftools norm -m -any "$out1" -Oz -o "$out2"

Get the sizes

    bcftools query -f '%CHROM\t%POS\t%TYPE\t%REF\t%ALT\n' "$out2" | awk 'BEGIN {OFS="\t"} {print $1, $2, $3, length($5) - length($4)}' > sv_sizes.txt

    bcftools view -i 'strlen(REF) < strlen(ALT)' "$out2" -Oz -o SV_insertion.vcf.gz
    bcftools view -i 'strlen(REF) > strlen(ALT)' "$out2" -Oz -o SV_deletion.vcf.gz

Then convert the vcf files into a fasta file using `scripts/convert_vcf_to_fasta.py`

Now that there's a fasta file, Repeatmasker can be run using:

    #!/bin/bash
    set -euo pipefail

    file="buffalo_pangenome_v3/pggb_concat_VCF/SV_location/SV_allindels.fa"
    lib="tools/repeatmasker/RepeatMasker/custom_lib/RepeatMaskerLib.h5.bostaurus-rm.withsats.fa"

    /hpcfs/groups/phoenix-hpc-avsci/Paulene_Pineda/tools/repeatmasker/RepeatMasker/RepeatMasker -pa 24 -gff -lib "$lib" -dir SV_indel "$file" -e ncbi
    echo "Done! :)"
