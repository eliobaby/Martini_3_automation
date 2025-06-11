#!/usr/bin/env python3
import pandas as pd
# 1) Read in both CSVs
# Adjust these filenames if yours differ exactly
df1 = pd.read_csv("Table_S1_from_653.csv", header=0)
df2 = pd.read_csv("logP_validation(deltaG_Working_10Jun25b).csv", header=0)

# 2) Identify columns by index positions (zero‐based)
# For df1 (Table_S1_from_653):
#   index 2 → CAS
#   index 3 → delta_g_exp
#   index 4 → delta_g_tristan
cas_col_df1           = df1.columns[2]
delta_g_exp_col       = df1.columns[3]
delta_g_tristan_col   = df1.columns[4]

# For df2 (logP_validation(653_Tristan)):
#   index 1 → CAS (compound)
#   index 4 → deltaG_Minh
cas_col_df2       = df2.columns[1]
delta_g_minh_col  = df2.columns[4]

# 3) Filter df2 to keep only rows where deltaG_Minh != 0
filtered_df2 = df2[df2[delta_g_minh_col] != 0]

# 4) Merge on CAS between filtered_df2 and df1
merged = pd.merge(
    filtered_df2[[cas_col_df2, delta_g_minh_col]],
    df1[[cas_col_df1, delta_g_exp_col, delta_g_tristan_col]],
    left_on=cas_col_df2,
    right_on=cas_col_df1,
    how="inner"
)

# 5) Select and rename columns for output
output = merged[[cas_col_df2, delta_g_minh_col, delta_g_exp_col, delta_g_tristan_col]].copy()
output.columns = ["CAS", "delta_G_Minh", "delta_G_exp", "delta_G_tristan"]

# 6) Write result to new CSV
output.to_csv("delta_testing_653.csv", index=False)
print(f"Wrote {len(output)} rows to delta_testing_653.csv")
