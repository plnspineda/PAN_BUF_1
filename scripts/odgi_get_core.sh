#!/bin/bash
touch core_accessory_unique.txt
for og in chr*_s100k/*.final.og; do
  base=$(basename "$og" .fa.f3ab435.11fba48.e859484.smooth.final.og)
  echo "$base" >> core_accessory_unique.txt
  odgi stats -i "$og" -a '#',0 >> core_accessory_unique.txt
done