#!/usr/bin/env python3

import os
import pandas as pd

# 1) Read in both CSVs
df1 = pd.read_csv("Table_S1_from_653.csv", header=0)
df2 = pd.read_csv("logP_validation(deltaG_all_11Jun25).csv", header=0)

# 2) Identify columns by index positions (zero‐based)
cas_col_df1         = df1.columns[2]  # CAS in df1
delta_g_exp_col     = df1.columns[3]
delta_g_tristan_col = df1.columns[4]

cas_col_df2       = df2.columns[1]  # CAS in df2
delta_g_minh_col  = df2.columns[4]

# 3) Filter df2 to keep only rows where deltaG_Minh != 0
filtered_df2 = df2[df2[delta_g_minh_col] != 0].copy()

# 4) Helpers to check folder existence and "S,S" exclusion
def has_working_dir(cas):
    return os.path.isdir(os.path.join("Half_Working", str(cas)))

def contains_ss(cas):
    """
    Returns True if Working/<cas>/<cas>.txt exists and
    any line from index 6 onward contains the substring "S,S".
    """
    txt_path = os.path.join("Working", str(cas), f"{cas}.txt")
    if not os.path.isfile(txt_path):
        return False
    with open(txt_path, 'r') as f:
        lines = f.readlines()
    for line in lines[6:]:
        if "S,S" in line:
            return True
    return False

# 5) Report and drop any CAS without a Working folder
missing_dirs = [
    cas for cas in filtered_df2[cas_col_df2].unique()
    if not has_working_dir(cas)
]
if missing_dirs:
    print("Warning: no Working/ folders found for these CAS values, skipping:")
    for cas in missing_dirs:
        print("  ", cas)
filtered_df2 = filtered_df2[filtered_df2[cas_col_df2].apply(has_working_dir)]
'''
# 6) Exclude any CAS whose .txt file contains "S,S" from line 6 onward
to_exclude_ss = [
    cas for cas in filtered_df2[cas_col_df2].unique()
    if contains_ss(cas)
]
if to_exclude_ss:
    print("\nExcluding CAS with 'S,S' in their .txt file:")
    for cas in to_exclude_ss:
        print("  ", cas)
filtered_df2 = filtered_df2[~filtered_df2[cas_col_df2].isin(to_exclude_ss)]
'''
# 7) Merge on CAS between filtered_df2 and df1
merged = pd.merge(
    filtered_df2[[cas_col_df2, delta_g_minh_col]],
    df1[[cas_col_df1, delta_g_exp_col, delta_g_tristan_col]],
    left_on=cas_col_df2,
    right_on=cas_col_df1,
    how="inner"
)

# 8) Select and rename columns for output
output = merged[[cas_col_df2, delta_g_minh_col, delta_g_exp_col, delta_g_tristan_col]].copy()
output.columns = ["CAS", "delta_G_Minh", "delta_G_exp", "delta_G_tristan"]

# 9) Write result to new CSV
output.to_csv("delta_testing_653.csv", index=False)
print(f"\nWrote {len(output)} rows to delta_testing_653.csv")
