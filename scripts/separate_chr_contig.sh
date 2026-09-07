#!/bin/bash
files=(*.fna)
for i in "${files[@]}"; do
    faidx -x "$i"
done