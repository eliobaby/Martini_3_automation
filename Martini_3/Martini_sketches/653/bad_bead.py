#!/usr/bin/env python3

import os
import csv
import pandas as pd

def count_beads_in_file(cas, counter_dict):
    """
    Given a CAS string and a dictionary (counter_dict),
    attempt to open Working/<cas>/<cas>.text; if not found, try Half_working/<cas>/<cas>.text.
    Parse lines from index 6 onward, count occurrences of the bead names (second token on each line),
    stopping when encountering a blank line. Update counter_dict in place.
    """
    # First attempt in Working/
    txt_path = os.path.join("Working", cas, f"{cas}.txt")
    if not os.path.isfile(txt_path):
        # Fallback to Half_working/
        txt_path = os.path.join("Half_working", cas, f"{cas}.txt")
        if not os.path.isfile(txt_path):
            # Neither path exists, skip this CAS
            return

    with open(txt_path, 'r') as f:
        lines = f.readlines()

    # Iterate from line index 6 onward
    for line in lines[6:]:
        if line.strip() == "":
            # blank or whitespace-only line → stop processing
            break
        parts = line.split()
        if len(parts) > 1:
            bead = parts[1]
            counter_dict[bead] = counter_dict.get(bead, 0) + 1


def main():
    # 1) Initialize empty dictionaries
    total = {}
    error = {}

    # 2) Read the CSV
    df = pd.read_csv("delta_testing_653.csv", header=0)

    # 3) Build the “total” dictionary:
    #    For each unique CAS in column index 0, count beads.
    unique_cas_all = df.iloc[:, 0].unique()
    for cas in unique_cas_all:
        count_beads_in_file(cas, total)

    # 4) Build the “error” dictionary:
    #    Compare column index 1 vs column index 2; if |diff| > 2, count beads.
    col_min = df.columns[1]
    col_exp = df.columns[2]

    mask = (df[col_min] - df[col_exp]).abs() > 1
    unique_cas_error = df.loc[mask, df.columns[0]].unique()

    for cas in unique_cas_error:
        count_beads_in_file(cas, error)

    # 5) Create bad_bead.csv:
    #    For each bead in total, get total_count, error_count (0 if missing), percent error = error/total.
    output_rows = []
    for bead, tot_count in total.items():
        err_count = error.get(bead, 0)
        percent_err = err_count / tot_count if tot_count != 0 else 0.0
        output_rows.append({
            "bead": bead,
            "error": err_count,
            "total": tot_count,
            "%error": percent_err
        })

    # 6) Sort by %error descending (highest error rate first)
    output_rows.sort(key=lambda r: r["%error"], reverse=True)

    # 7) Write to bad_bead.csv
    fieldnames = ["bead", "error", "total", "%error"]
    with open("bad_bead_2.csv", "w", newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in output_rows:
            writer.writerow(row)

    print(f"Wrote {len(output_rows)} beads to bad_bead_2.csv")


if __name__ == "__main__":
    main()
