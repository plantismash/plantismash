import json

# Define your target gene IDs (strip "_logo.png")
motif_ids = [
    "AT1G28300",
    "AT1G63030",
    "AT2G28550",
    "AT2G46510",
    "AT3G14180"
]

# Load the full motif JSON
with open("Athaliana.json", "r") as f:
    full_data = json.load(f)

# Filter only selected motifs
filtered_data = {key: full_data[key] for key in motif_ids if key in full_data}

# Save the new filtered JSON
with open("Athaliana_subset.json", "w") as f:
    json.dump(filtered_data, f, indent=4)

print(f"✅ Saved {len(filtered_data)} motifs to Athaliana_subset.json")
