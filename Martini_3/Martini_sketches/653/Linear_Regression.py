#!/usr/bin/env python3
import os
import pandas as pd
import numpy as np

def get_beads(cas):
    """Return list of beads for a given CAS."""
    for root in ("Working","Half_working"):
        for ext in (".text",".txt"):
            path = os.path.join(root, cas, f"{cas}{ext}")
            if os.path.isfile(path):
                with open(path) as f:
                    lines = f.readlines()[6:]
                beads = []
                for L in lines:
                    if not L.strip(): break
                    parts = L.split()
                    if len(parts)>1:
                        beads.append(parts[1])
                return beads
    return []

def main():
    # 1) Load delta_testing_653.csv
    df = pd.read_csv("delta_testing_653.csv", header=0)
    # compute per‐molecule absolute error
    df["abs_error"] = (df.iloc[:,1] - df.iloc[:,2]).abs()

    # 2) Gather bead lists & vocabulary
    cas_list = df.iloc[:,0].astype(str).tolist()
    cas_to_beads = {}
    vocab = set()
    for cas in cas_list:
        beads = get_beads(cas)
        cas_to_beads[cas] = beads
        vocab.update(beads)
    vocab = sorted(vocab)

    # 3) Build design matrix C (M×N) and target E (M)
    M = len(cas_list)
    N = len(vocab)
    C = np.zeros((M, N), dtype=float)
    for i, cas in enumerate(cas_list):
        counts = pd.value_counts(cas_to_beads[cas])
        for j, bead in enumerate(vocab):
            C[i, j] = counts.get(bead, 0)
    E = df["abs_error"].values

    # 4) Solve least‐squares: x = argmin ||C x − E||
    x, residuals, rank, s = np.linalg.lstsq(C, E, rcond=None)

    # 5) Save bead contributions
    out = pd.DataFrame({
        "bead": vocab,
        "estimated_error": x
    }).sort_values("estimated_error", ascending=False)
    out.to_csv("bead_error_contributions.csv", index=False)
    print(f"Wrote {len(out)} beads to bead_error_contributions.csv")

if __name__=="__main__":
    main()
