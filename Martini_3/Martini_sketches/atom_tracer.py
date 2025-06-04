from rdkit import Chem
from rdkit.Chem import rdDepictor
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Geometry import Point3D
from PIL import Image
import io

# 1) Read your SMILES, build the mol, drop explicit H's, compute 2D coords
smiles = input("Smiles: ").strip()
mol = Chem.MolFromSmiles(smiles)
mol = Chem.RemoveHs(mol)
rdDepictor.Compute2DCoords(mol)

# 2) Scale the drawing outwards so bonds look longer
conf = mol.GetConformer()
# %%
scale = 2  # e.g. 150% size

for i in range(mol.GetNumAtoms()):
    p = conf.GetAtomPosition(i)                       # this is already a Point3D(x,y,z)
    conf.SetAtomPosition(i,
                         Point3D(p.x * scale,
                                 p.y * scale,
                                 p.z))                # keep z the same

# 3) Set up a large drawing surface
width, height = 2000, 2000
drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
opts = drawer.drawOptions()
opts.annotationFontScale = 0.5
opts.fontSize          = 8

# 4) Override each atom label with “Element+Index”
for atom in mol.GetAtoms():
    idx    = atom.GetIdx()
    symbol = atom.GetSymbol()
    opts.atomLabels[idx] = f"{symbol}{idx}"

# 5) Draw & finish
drawer.DrawMolecule(mol)
drawer.FinishDrawing()

# 6) Convert to PIL Image and show
png = drawer.GetDrawingText()
img = Image.open(io.BytesIO(png))
img.show()
