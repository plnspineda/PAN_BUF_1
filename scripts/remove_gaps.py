import sys
import os
import subprocess

input_fasta = sys.argv[1]
output_ungapped = os.path.splitext(input_fasta)[0] + "_ungapped.fa"
mapping_file = os.path.splitext(input_fasta)[0] + "_mapping.txt"
tmp_file = "tmp.fa"

with open(tmp_file, "w") as outfile:
    subprocess.run(["seqtk", "cutN", "-n", "3", input_fasta], stdout=outfile)

with open(tmp_file, "r") as infile, open(output_ungapped, "w") as outfile, open(mapping_file, "w") as mapfile:
    index = 0
    for line in infile:
        if line.startswith(">"):
            parts = line.strip().split(":")
            header = parts[0]
            coordinates = parts[1]
            new_header = f"{header}#piece{index}"
            outfile.write(f"{new_header}\n")
            mapfile.write(f"{index} {coordinates}\n")
            index += 1
        else:
            outfile.write(line)

os.remove(tmp_file)