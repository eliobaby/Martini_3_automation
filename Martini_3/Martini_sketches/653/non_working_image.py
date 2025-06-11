#!/usr/bin/env python3
import os
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw
from PIL import Image

def main():
    # 1. Read the list of non-working CAS identifiers
    with open("non-working.txt") as f:
        cas_list = [line.strip() for line in f if line.strip()]

    # 2. Load the CSV with columns: 0='#', 1='SMILES', 2='CAS', ...
    df = pd.read_csv("Table_S1_from_653.csv", dtype=str)
    cas_col    = df.columns[2]
    smiles_col = df.columns[1]

    # 3. Collect RDKit Mol objects for each CAS
    mols = []
    for cas in cas_list:
        match = df[df[cas_col] == cas]
        if match.empty:
            print(f"Warning: CAS {cas} not found in the CSV.")
            continue
        smiles = match[smiles_col].iloc[0]
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            print(f"Warning: invalid SMILES for CAS {cas}: {smiles}")
            continue
        mols.append(mol)

    # 4. Define grid layout and image parameters
    cols, rows = 6, 10
    dpi_x, dpi_y = 50, 50
    # each sub-image is 2"×2", so pixel size = inches × DPI
    sub_w = int(2 * dpi_x)
    sub_h = int(2 * dpi_y)

    # 5. Create a blank white canvas
    mosaic_w = cols * sub_w
    mosaic_h = rows * sub_h
    mosaic = Image.new("RGB", (mosaic_w, mosaic_h), "white")

    # 6. Draw and paste each molecule into the grid
    for idx, mol in enumerate(mols):
        if idx >= cols * rows:
            break
        img = Draw.MolToImage(mol, size=(sub_w, sub_h))
        col = idx % cols
        row = idx // cols
        mosaic.paste(img, (col * sub_w, row * sub_h))

    # 7. Save the final mosaic
    out_name = "non_working_molecules.png"
    mosaic.save(out_name, dpi=(dpi_x, dpi_y))
    print(f"Saved {len(mols)} structures to {out_name} at grid {cols}×{rows}, DPI={dpi_x}×{dpi_y}")

if __name__ == "__main__":
    main()
