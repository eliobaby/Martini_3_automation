from rdkit import Chem
from setup import get_atom_properties, connectivity_matrix, output_comb_polymer_visual
from setup_canonical import remove_special_characters, transform_smiles
from matrix_maker import comb_check, matrix_size_calc, backbone_maker, sidechain_maker, bead_empty_matrix
from test1_C_only import get_c_only_dict
from algorithm import triple_bond_assignment, full_bead_assignment

#Input SMILES string
smiles = input("Smiles: ")

mol = Chem.MolFromSmiles(smiles)

atom_properties = get_atom_properties(mol)

connectivity_matrix = connectivity_matrix(mol, len(atom_properties))

compound_name = input("Name: ")

smiles = remove_special_characters(smiles)
output_comb_polymer_visual(smiles, compound_name)

bead_dict = get_c_only_dict()

#Check if the Smiles string is valid
smiles = transform_smiles(smiles)
comb_check(smiles)

#Create the matrix needed for the actual algorithm
matrix_empty = matrix_size_calc(smiles)
matrix_backbone = backbone_maker(smiles, matrix_empty)
matrix_full = sidechain_maker(smiles, matrix_backbone)

#Create an empty matrix to store the results
bead_empty = bead_empty_matrix(matrix_empty)

matrix_full, bead_empty = full_bead_assignment(matrix_full, bead_empty)
matrix_full, bead_empty = triple_bond_assignment(matrix_full, bead_empty)

print(bead_empty)