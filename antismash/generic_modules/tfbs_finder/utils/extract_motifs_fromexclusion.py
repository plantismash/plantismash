#!/usr/bin/env python3
import json
from pathlib import Path

SOURCE_JSON   = Path("../data/Athaliana_motifs.filtered.json")
EXCLUDE_LIST  = Path("../data/motifs_excluded.txt")
OUT_JSON      = Path("../data/Athaliana_motifs.manualfromexcluded.json")

# 1) read motif IDs to exclude
to_exclude = []
for line in EXCLUDE_LIST.read_text().splitlines():
    line = line.strip()
    if not line or line.startswith("#"):
        continue
    to_exclude.append(line)

exclude_set = set(to_exclude)

# 2) load source JSON (expects a dict keyed by motif IDs)
data = json.loads(SOURCE_JSON.read_text())

# 3) filter: keep only motifs NOT in exclude_set
subset = {k: v for k, v in data.items() if k not in exclude_set}

# 4) write output (pretty + stable key order)
OUT_JSON.write_text(json.dumps(subset, indent=2, sort_keys=True))

# 5) simple report
n_total   = len(data)
n_removed = n_total - len(subset)
n_kept    = len(subset)

missing = [m for m in to_exclude if m not in data]

print(f"Total motifs in source JSON: {n_total}")
print(f"Excluded motifs: {n_removed}")
print(f"Remaining motifs: {n_kept} → {OUT_JSON}")
if missing:
    print("IDs in exclude list not found in source JSON:", ", ".join(missing))