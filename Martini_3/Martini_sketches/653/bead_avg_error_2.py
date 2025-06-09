#!/usr/bin/env python3

import os
import csv
import pandas as pd

def update_bead_errors(cas, error_value, bead_dict):
    """
    Given a CAS string, an absolute error value, and the bead_dict,
    open Working/<cas>/<cas>.text (or .txt); if not found, try Half_working/<cas>/<cas>.text (or .txt).
    From line index 6 onward, read each line until a blank line; split by whitespace and take parts[1] as bead name.
    For each bead, update bead_dict[bead] = [running_avg_error, count, set_of_cas]:
      - If bead not in bead_dict, initialize [error_value, 1, {cas}]
      - Else, let [prev_avg, prev_count, prev_set] = bead_dict[bead].
           new_count = prev_count + 1
           new_avg = (prev_avg * prev_count + error_value) / new_count
           new_set = prev_set ∪ {cas}
           store [new_avg, new_count, new_set]
    """
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
        return

    with open(found_path, "r") as f:
        lines = f.readlines()

    for line in lines[6:]:
        if not line.strip():
            break
        parts = line.split()
        if len(parts) > 1:
            bead = parts[1]
            if bead not in bead_dict:
                bead_dict[bead] = [error_value, 1, {cas}]
            else:
                prev_avg, prev_count, prev_set = bead_dict[bead]
                new_count = prev_count + 1
                new_avg = (prev_avg * prev_count + error_value) / new_count
                prev_set.add(cas)
                bead_dict[bead] = [new_avg, new_count, prev_set]


def main():
    bead_dict = {}

    csv_path = "delta_testing_653.csv"
    if not os.path.isfile(csv_path):
        print(f"Error: {csv_path} not found in the current directory. Exiting.")
        return

    df = pd.read_csv(csv_path, header=0)

    for _, row in df.iterrows():
        cas = str(row.iloc[0])
        delta_minh = float(row.iloc[1])
        delta_exp = float(row.iloc[2])
        abs_error = (delta_minh - delta_exp)
        update_bead_errors(cas, abs_error, bead_dict)

    # Prepare output rows
    output_rows = []
    for bead, (avg_err, _, cas_set) in bead_dict.items():
        cas_list = ";".join(sorted(cas_set))
        output_rows.append({
            "bead": bead,
            "avg_error": avg_err,
            "CAS_list": cas_list
        })

    # Sort by avg_error descending
    output_rows.sort(key=lambda r: r["avg_error"], reverse=True)

    # Write to bead_avg_error.csv
    out_csv = "bead_avg_error_2.csv"
    fieldnames = ["bead", "avg_error", "CAS_list"]
    with open(out_csv, "w", newline="") as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in output_rows:
            writer.writerow(row)

    print(f"Wrote {len(output_rows)} beads (with CAS lists) to {out_csv}")


if __name__ == "__main__":
    main()
