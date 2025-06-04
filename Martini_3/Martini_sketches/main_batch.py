import copy
import sys
from rdkit import Chem
from setup import get_atom_properties, connectivity_matrix, duplicate_ring_number_handler
from martini_3_dictionary import get_m3_dict
from mapping_scheme import map_molecule, parse_smiles
from algorithm import map_martini_beads
from text_file import gro_maker, itp_maker, export_bead_mapping

def process_fragment(compound_name: str, smiles: str, suffix: str = ""):
    """
    Process a single SMILES fragment and write its bead mapping.
    """
    mol = Chem.MolFromSmiles(smiles)
    atom_props = get_atom_properties(mol)
    final = ["" for _ in atom_props]
    conn_mat = connectivity_matrix(mol, len(atom_props))
    dict3 = get_m3_dict()
    smiles_token = parse_smiles(smiles)
    mapping = map_molecule(smiles_token, conn_mat, atom_props)
    mapping_copy = copy.deepcopy(mapping)
    final = map_martini_beads(mapping, final, dict3)
    name = compound_name + suffix
    beads_data, connections_data = export_bead_mapping(final, mapping_copy, smiles, name)
    gro_maker(beads_data, name)
    itp_maker(beads_data, connections_data, name)


def main():
    compound_name = input().strip()
    smiles        = input().strip()
    smiles        = duplicate_ring_number_handler(smiles)

    # If there are multiple fragments, pick only the longest one
    if '.' in smiles:
        frags   = smiles.split('.')
        longest = max(frags, key=len)
        try:
            process_fragment(compound_name, longest)
        except Exception as e:
            print(f"Error: failed on largest fragment '{longest}': {e}")
            sys.exit(1)
        sys.exit(0)

    # Single‐fragment case
    try:
        process_fragment(compound_name, smiles)
    except Exception as e:
        print(f"Error: {e}")
        sys.exit(1)
        
if __name__ == '__main__':
    main()
