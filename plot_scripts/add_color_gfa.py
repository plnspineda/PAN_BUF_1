#!/usr/bin/env python3

import sys
import csv
from collections import defaultdict

if len(sys.argv) != 3:
    print(f"Usage: {sys.argv[0]} graph.gfa bandage_colors.csv")
    sys.exit(1)

gfa = sys.argv[1]
outfile = sys.argv[2]

# node -> set(path names)
node_paths = defaultdict(set)

# node -> depth
depth = {}

with open(gfa) as f:
    for line in f:
        if line.startswith("S"):
            fields = line.rstrip().split("\t")
            node = fields[1]

            d = None
            for tag in fields[3:]:
                if tag.startswith("DP:f:"):
                    d = float(tag.split(":")[-1])
                elif tag.startswith("dp:i:"):
                    d = int(tag.split(":")[-1])
                elif tag.startswith("DP:i:"):
                    d = int(tag.split(":")[-1])

            if d is not None:
                depth[node] = d

        elif line.startswith("P"):
            fields = line.rstrip().split("\t")
            pname = fields[1]

            nodes = fields[2].split(",")

            for n in nodes:
                node = n.rstrip("+-")
                node_paths[node].add(pname)

with open(outfile, "w", newline="") as out:
    writer = csv.writer(out)
    writer.writerow(["Name", "Color", "Label"])

    for node in sorted(node_paths.keys(), key=lambda x: int(x) if x.isdigit() else x):

        # Priority 1: depth == 9
        if depth.get(node) == 9:
            color = "#00AA00"     # green

        else:
            paths = node_paths[node]

            has_sw = any(p.startswith("SW") for p in paths)
            has_rv = any(
                p.startswith("RV") or p.startswith("ARV")
                for p in paths
            )

            if has_sw and has_rv:
                color = "#8A2BE2"     # violet
            elif has_sw:
                color = "#FF69B4"     # pink
            elif has_rv:
                color = "#1E90FF"     # blue
            else:
                color = "#B0B0B0"     # gray

        label = ",".join(sorted(node_paths[node]))
        writer.writerow([node, color, label])

print(f"Wrote {outfile}")
