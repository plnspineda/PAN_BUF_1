#!/usr/bin/env python3

import sys
import gzip
import argparse
import re

def process_vcf(vcf_file):
    if vcf_file.endswith(".gz"):
        with gzip.open(vcf_file, "rt") as f:
            for line in f:
                process_line(line)
    else:
        with open(vcf_file, "r") as f:
            for line in f:
                process_line(line)

def process_line(line):
    if line.startswith("##") or line.startswith("#CHROM"):
        return

    cols = line.strip().split("\t")
    if len(cols) < 5:
        return

    chrom = cols[0]
    pos = cols[1]
    variant_id = cols[2]
    ref_seq = cols[3]
    alt_seqs = cols[4].split(",")

    chrom_number = chrom
    start = pos
    end = pos
    suffix = "0"
    
    match = re.match(r"(\d+)_(\d+)_(\d+)_([^\s_]+)", variant_id)
    if match:
        chrom_number = match.group(1)
        start = match.group(2)
        end = match.group(3)
        suffix = match.group(4)

    ref_length = len(ref_seq)
    alt_lengths = [len(alt) for alt in alt_seqs]
    max_length = max([ref_length] + alt_lengths)

    variant_type = "ins" if ref_length == 1 else "del"

    if max_length == ref_length:
        seq_type = "REF"
        longest_sequence = ref_seq
    else:
        seq_type = "ALT"
        longest_sequence = alt_seqs[alt_lengths.index(max_length)]

    fasta_id = f"{variant_type}_{chrom_number}_{start}_{suffix}_{max_length}"
    print(f">{fasta_id}")
    print(longest_sequence)

def main():
    parser = argparse.ArgumentParser(
        description="Extract the longest sequence (REF or ALT) for each variant in a VCF file and output as FASTA."
    )
    parser.add_argument("-v", "--vcf", required=True, help="Input VCF file (.vcf or .vcf.gz)")

    args = parser.parse_args()
    process_vcf(args.vcf)

if __name__ == "__main__":
    main()
