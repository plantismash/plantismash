#!/usr/bin/env python3
import json
from pathlib import Path

SOURCE_JSON = Path("../data/Athaliana_motifs.filtered.json")
MOTIF_LIST  = Path("motifs_manuallycurated_names.txt")
OUT_JSON    = Path("../data/Athaliana_motifs.subset.json")

# 1) read motif IDs
motifs = []
for line in MOTIF_LIST.read_text().splitlines():
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    motifs.append(line)

motif_set = set(motifs)

# 2) load source JSON (expects a dict keyed by motif IDs)
data = json.loads(SOURCE_JSON.read_text())

# 3) filter
subset = {k: v for k, v in data.items() if k in motif_set}

# 4) write output (pretty + stable key order)
OUT_JSON.write_text(json.dumps(subset, indent=2, sort_keys=True))

# 5) simple report
missing = [m for m in motifs if m not in data]
print(f"Kept {len(subset)} / {len(motifs)} motifs → {OUT_JSON}")
if missing:
    print("Not found in source JSON:", ", ".join(missing))