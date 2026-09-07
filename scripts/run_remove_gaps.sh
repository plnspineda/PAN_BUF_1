#!/bin/bash

module use /apps/modules/all
module load Python/3.11.3-GCCcore-12.3.0
module load seqtk/1.3-GCC-11.2.0

for i in *fna; do
    echo "python remove_gaps.py $i"
    python remove_gaps.py "$i"
done