import re
import pandas as pd
from pypdf import PdfReader

# 1) Read the PDF directly by filename
reader = PdfReader("653.pdf")

# 2) Concatenate all text
all_text = []
for page in reader.pages:
    text = page.extract_text()
    if text:
        all_text.append(text)
full_text = "\n".join(all_text)

# 3) Regex for Table S1:
#    - entry number (integer)
#    - SMILES (no whitespace)
#    - CAS  (pattern ddd-ddd-d)
#    - delta G experimental (float, can be negative)
#    - delta G “Tristan” (float, can be negative)
pattern = re.compile(
    r'^\s*(\d+)\s+'                              # 1) entry number
    r'(\S+)\s+'                                  # 2) SMILES
    r'(\d{1,7}-\d{1,7}-\d+)\s+'                  # 3) CAS
    r'(-?\d+\.\d+)\s+'                           # 4) ∆Gexp (signed float)
    r'(-?\d+\.\d+)'                              # 5) ∆Gtristan (signed float)
    , re.MULTILINE
)

rows = pattern.findall(full_text)

# 4) Build DataFrame, rename columns, sort by CAS, and write CSV
df = pd.DataFrame(rows, columns=['#','SMILES','CAS','delta_g_exp','delta_g_tristan'])

# Convert numeric columns to float
df['delta_g_exp'] = df['delta_g_exp'].astype(float)
df['delta_g_tristan'] = df['delta_g_tristan'].astype(float)

# Sort by the CAS column (lexicographically)
df = df.sort_values(by='CAS').reset_index(drop=True)

# Write out
df.to_csv("Table_S1_from_653.csv", index=False)
print(f"Wrote {len(df)} rows (including ΔG columns) to Table_S1_from_653.csv")
