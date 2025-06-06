#!/usr/bin/env python3

import os
import csv
import pandas as pd
import numpy as np

def update_bead_errors(cas, error_value, bead_dict):
    """
    Given a CAS string, an absolute error value, and the bead_dict,
    open Working/<cas>/<cas>.text (or .txt); if not found, try Half_working/<cas>/<cas>.text (or .txt).
    From line index 6 onward, read each line until a blank line; split by whitespace and take parts[1] as bead name.
    For each bead, update bead_dict[bead] = [running_avg_error, count]:
      - If bead not in bead_dict, initialize [error_value, 1]
      - Else, let [prev_avg, prev_count] = bead_dict[bead].
           new_count = prev_count + 1
           new_avg = (prev_avg * prev_count + error_value) / new_count
           store [new_avg, new_count]
    """
    # Try both Working and Half_working, and both extensions (.text, .txt)
    roots = ["Working", "Half_working"]
    exts = [".text", ".txt"]
    found_path = None

    for root in roots:
        for ext in exts:
            candidate = os.path.join(root, cas, f"{cas}{ext}")
            if os.path.isfile(candidate):
                found_path = candidate
                break
        if found_path:
            break

    if not found_path:
        # If no file was found under Working or Half_working, skip
        return

    with open(found_path, "r") as f:
        lines = f.readlines()

    # Iterate from line index 6 onward
    for line in lines[6:]:
        if line.strip() == "":
            # blank/whitespace-only → stop reading beads
            break
        parts = line.split()
        if len(parts) > 1:
            bead = parts[1]
            if bead not in bead_dict:
                # Initialize with [average_error, count]
                bead_dict[bead] = [error_value, 1]
            else:
                prev_avg, prev_count = bead_dict[bead]
                new_count = prev_count + 1
                new_avg = (prev_avg * prev_count + error_value) / new_count
                bead_dict[bead] = [new_avg, new_count]


def main():
    # 1) Initialize empty bead dictionary: { bead_name: [avg_error, count] }
    bead_dict = {}

    # 2) Read delta_testing_653.csv
    csv_path = "delta_testing_653.csv"
    if not os.path.isfile(csv_path):
        print(f"Error: {csv_path} not found in the current directory. Exiting.")
        return

    df = pd.read_csv(csv_path, header=0)

    # 3) Iterate over each row exactly once
    #    Columns: index 0 = CAS, index 1 = delta_G_Minh, index 2 = delta_G_exp
    for idx, row in df.iterrows():
        cas = str(row.iloc[0])
        delta_minh = float(row.iloc[1])
        delta_exp = float(row.iloc[2])
        abs_error = abs(delta_minh - delta_exp)

        # Update bead_dict for this CAS using abs_error
        update_bead_errors(cas, abs_error, bead_dict)

    # 4) After processing all rows, create bead_avg_error.csv
    output_rows = []
    for bead, (avg_err, count) in bead_dict.items():
        output_rows.append({
            "bead": bead,
            "avg_error": avg_err
        })

    # 5) Sort by avg_error descending (highest average error first)
    output_rows.sort(key=lambda r: r["avg_error"], reverse=True)

    # 6) Write to bead_avg_error.csv
    out_csv = "bead_avg_error.csv"
    fieldnames = ["bead", "avg_error"]
    with open(out_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in output_rows:
            writer.writerow(row)

    print(f"Wrote {len(output_rows)} beads to {out_csv}")


if __name__ == "__main__":
    main()
