# Benchmarking of SVs from pangenome graph

Use `scripts/control_SVIM_array` to run benchmarking with SVIM and PGGB variants.

These are the parameters used for SURVIVOR and truvari

    SURVIVOR merge vcf_list.txt 500 1 1  0  50  2 \
    merged_SVs.buffalo_pop.noINVBND.vcf

    truvari bench \
    -b pggb.norm.vcf.gz \
    -c merged_SVs.buffalo_pop.INSDEL.sorted.vcf.gz \
    -o truvari_pggb_base_unfiltered \
    --refdist 500 \
    --pctsize 0.5 \
    --pctseq 0.0 \
    --sizemin 50 \
    --passonly

