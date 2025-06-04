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

# 3) Regex for Table S1: entry number, SMILES, CAS
pattern = re.compile(
    r'^\s*(\d+)\s+'            # entry number
    r'(\S+)\s+'                # SMILES (no spaces)
    r'(\d{1,7}-\d{1,7}-\d+)'   # CAS number
    , re.MULTILINE
)
rows = pattern.findall(full_text)

# 4) Build DataFrame and write CSV
df = pd.DataFrame(rows, columns=['#','SMILES','CAS'])
df.to_csv("Table_S1_from_653.csv", index=False)
print(f"Wrote {len(df)} rows to Table_S1_from_653.csv")
