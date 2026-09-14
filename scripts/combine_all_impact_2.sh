awk -F'\t' -v OFS='\t' '
    /^#/ { next } 
    {
        if (!seen[$2]++) {
            split($14, a, ";");
            split(a[1], b, "=");

            label = FILENAME;
            sub(/.*_r1M\./, "", label);
            sub(/\.txt$/, "", label);
            gsub(/\./, "_", label);

            # Added $7 right here
            print $2, b[2], $7, label
        }
    }
' ARVUO_pggb_all_chr_p95_s10k_chr8_p93_vcfbub_L0_r1M.*.txt > impact_all_2.txt
