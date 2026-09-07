#!/bin/bash

out="pangenome_size_stats.txt"
echo -e "ID\t#length\t#nodes\t#edges\t#paths\t#steps\t#segments" > "$out"

for og in chr*s10k/*smooth.final.og; do
  base=$(basename "$og")
  base=${base%%.*}
  stats=$(odgi stats -i "$og" -S | tail -n 1)

  echo -e "${base}\t${stats}" >> "$out"
done
