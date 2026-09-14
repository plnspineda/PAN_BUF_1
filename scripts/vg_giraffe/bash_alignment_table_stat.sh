#!/bin/bash

# Output file
in="$1"
base="$(basename $in)"
out="alignment_summary_$base.tsv"

# Print header
echo -e "File\tTotal_alignments\tTotal_primary\tTotal_secondary\tTotal_aligned\tTotal_perfect\tTotal_gapless\tTotal_paired\tTotal_properly_paired\tAS_mean\tAS_median\tAS_stdev\tAS_max\tMQ_mean\tMQ_median\tMQ_stdev\tMQ_max\tInsertions_bp\tInsertion_events\tDeletions_bp\tDeletion_events\tSubstitutions_bp\tSubstitution_events\tSoftclips_bp\tSoftclip_events\tTotal_time_sec\tSpeed_reads_per_sec" > "$out"

# Loop through files (change pattern if needed)
find "$in" -type f -name "*stats" | while read -r f; do
    awk -v file="$f" '
    BEGIN {
        OFS="\t"
    }

    /Total alignments:/        {ta=$3}
    /Total primary:/           {tp=$3}
    /Total secondary:/         {ts=$3}
    /Total aligned:/           {tal=$3}
    /Total perfect:/           {tperf=$3}
    /Total gapless/            {tgap=$5}
    /Total paired:/            {tpaired=$3}
    /Total properly paired:/   {tproper=$4}

    /Alignment score:/ {
        as_mean=$4
        as_median=$6
        as_stdev=$8
        as_max=$10
        gsub(",", "", as_mean)
        gsub(",", "", as_median)
        gsub(",", "", as_stdev)
        gsub(",", "", as_max)
    }

    /Mapping quality:/ {
        mq_mean=$4
        mq_median=$6
        mq_stdev=$8
        mq_max=$10
        gsub(",", "", mq_mean)
        gsub(",", "", mq_median)
        gsub(",", "", mq_stdev)
        gsub(",", "", mq_max)
    }

    /Insertions:/      {ins_bp=$2; ins_ev=$5}
    /Deletions:/       {del_bp=$2; del_ev=$5}
    /Substitutions:/   {sub_bp=$2; sub_ev=$5}
    /Softclips:/       {soft_bp=$2; soft_ev=$5}

    /Total time:/      {time=$3}
    /Speed:/           {speed=$2}

    END {
        print file, ta, tp, ts, tal, tperf, tgap, tpaired, tproper,
              as_mean, as_median, as_stdev, as_max,
              mq_mean, mq_median, mq_stdev, mq_max,
              ins_bp, ins_ev,
              del_bp, del_ev,
              sub_bp, sub_ev,
              soft_bp, soft_ev,
              time, speed
    }
    ' "$f" >> "$out"

done

echo "Done. Output written to $out"
