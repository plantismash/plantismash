#!/usr/bin/env python3

import json

INPUT_FILE = "Athaliana.json"
OUTPUT_FILE = "Athaliana_transposed.json"

def transpose_pwm(pwm):
    """Transpose PWM from Nx4 to 4xN format."""
    return list(map(list, zip(*pwm)))

def main():
    with open(INPUT_FILE, "r") as infile:
        motifs = json.load(infile)

    for motif in motifs.values():
        pwm = motif.get("pwm")
        if pwm and len(pwm[0]) == 4:  # Only transpose if not already transposed
            motif["pwm"] = transpose_pwm(pwm)

    with open(OUTPUT_FILE, "w") as outfile:
        json.dump(motifs, outfile, indent=4)

    print(f"✅ Transposed {len(motifs)} PWMs and saved to '{OUTPUT_FILE}'.")

if __name__ == "__main__":
    main()
