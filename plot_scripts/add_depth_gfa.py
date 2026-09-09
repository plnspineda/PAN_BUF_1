#!/usr/bin/env python3

import sys
from collections import defaultdict

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} input.gfa output.gfa")
    sys.exit(1)

input_gfa = sys.argv[1]
output_gfa = sys.argv[2]

# node -> set of path names
node_paths = defaultdict(set)

# Read all paths
with open(input_gfa) as f:
    for line in f:
        if line.startswith("P\t"):
            fields = line.rstrip().split("\t")

            path_name = fields[1]
            nodes = fields[2].split(",")

            seen = set()
            for node in nodes:
                node = node.rstrip("+-")
                seen.add(node)

            for node in seen:
                node_paths[node].add(path_name)

# Rewrite GFA
with open(input_gfa) as fin, open(output_gfa, "w") as fout:
    for line in fin:
        if line.startswith("S\t"):
            fields = line.rstrip().split("\t")
            node = fields[1]

            # remove old tags if present
            fields = [x for x in fields
                      if not x.startswith("DP:i:")
                      and not x.startswith("HN:Z:")]

            paths = sorted(node_paths.get(node, []))
            depth = len(paths)

            fields.append(f"DP:i:{depth}")
            fields.append(f"HN:Z:{','.join(paths)}")

            fout.write("\t".join(fields) + "\n")
        else:
            fout.write(line)
