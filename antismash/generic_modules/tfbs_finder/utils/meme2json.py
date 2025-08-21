#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import json
import os

BASES = ['A', 'C', 'G', 'T']

def compute_min_max_scores(pwm):
    """Compute the min and max possible scores for a PWM."""
    min_score = sum(min(col) for col in pwm)
    max_score = sum(max(col) for col in pwm)
    return round(min_score, 4), round(max_score, 4)

def compute_consensus(pwm):
    """Compute the consensus sequence from the PWM."""
    consensus = ""
    for col in pwm:
        max_index = col.index(max(col))
        consensus += BASES[max_index]
    return consensus

def parse_meme_file(meme_file):
    motifs = {}
    current_motif = None
    pwm = []

    with open(meme_file, 'r') as f:
        for line in f:
            line = line.strip()

            if not line or line.startswith("MEME") or "ALPHABET=" in line or "strands:" in line:
                continue

            if line.startswith("MOTIF"):
                parts = line.split()
                if current_motif and pwm:
                    current_motif["pwm"] = pwm
                    current_motif["min_score"], current_motif["max_score"] = compute_min_max_scores(pwm)
                    current_motif["consensus"] = compute_consensus(pwm)
                    motifs[current_motif["name"]] = current_motif
                pwm = []
                if len(parts) >= 2:
                    current_motif = {
                        "name": parts[1],
                        "description": "Unknown transcription factor",
                        "species": "Arabidopsis thaliana",
                        "link": "",
                        "consensus": "",
                        "max_score": 0.0,
                        "min_score": 0.0,
                        "pwm": []
                    }

            elif current_motif and line and line[0].isdigit():
                try:
                    values = list(map(float, line.split()))
                    if len(values) == 4:
                        pwm.append(values)
                except ValueError:
                    continue

            elif line.startswith("URL") and current_motif:
                current_motif["link"] = line.split(" ", 1)[1]

    if current_motif and pwm:
        current_motif["pwm"] = pwm
        current_motif["min_score"], current_motif["max_score"] = compute_min_max_scores(pwm)
        current_motif["consensus"] = compute_consensus(pwm)
        motifs[current_motif["name"]] = current_motif

    return motifs

if __name__ == "__main__":
    input_meme_file = "Ath_TF_binding_motifs.meme"
    output_json_file = "Athaliana.json"

    motifs_dict = parse_meme_file(input_meme_file)

    with open(output_json_file, 'w') as json_file:
        json.dump(motifs_dict, json_file, indent=4)

    print(f"✅ JSON file '{output_json_file}' created with {len(motifs_dict)} motifs.")
