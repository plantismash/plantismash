#!/usr/bin/env python3
import os

# Directory where your PNG files are
directory = "../data/logo_included/"

# Output file
outfile = "motifs_manuallycurated_names.txt"

motifs = []
for filename in os.listdir(directory):
    if filename.endswith("_logo.png"):
        # remove the suffix
        motif = filename.replace("_logo.png", "")
        motifs.append(motif)

# Write to file
with open(outfile, "w") as f:
    for motif in sorted(motifs):
        f.write(motif + "\n")

print(f"Extracted {len(motifs)} motifs into {outfile}")