import numpy as np

def get_atom_properties(mol):
    #Create the array of arrays with atomic properties
    atom_properties = []
    for atom in mol.GetAtoms():
        mass = atom.GetMass()
        neighbor = atom.GetNeighbors()
        atom_data = [
            atom.GetSymbol(),            # Atom type (e.g., C, O, N)
            atom.GetMass(),              # Atomic mass
            atom.GetFormalCharge(),      # Formal charge
            atom.IsInRing(),             # Part of a ring?
            any(mass > 1.1 for neighbor in neighbor) and atom.GetDegree() == 1
            # Is it an edge node? (connects to only 1 other large node)
        ]
        atom_properties.append(atom_data)
    return atom_properties

def connectivity_matrix(mol, length):
    #Create the NxN connectivity matrix
    connectivity_matrix = np.zeros((length, length), dtype=int)

    for bond in mol.GetBonds():
        atom1_idx = bond.GetBeginAtomIdx()
        atom2_idx = bond.GetEndAtomIdx()
        bond_order = bond.GetBondTypeAsDouble()  # Can be 1 for single bond, 2 for double bond, etc.
        connectivity_matrix[atom1_idx][atom2_idx] = bond_order
        connectivity_matrix[atom2_idx][atom1_idx] = bond_order
    return connectivity_matrix

def output_comb_polymer_visual(smiles, compound_name):
    """
    Create a text file representing the comb polymer from its canonical SMILES.
    
    The SMILES is assumed to be written with:
      - The backbone outside parentheses.
      - One or more sidechains inside parentheses.
      - Bond symbols: "-" for single bonds (if missing, inserted automatically),
                      "=" for double bonds, and "#" for triple bonds.
    
    The output grid is built as follows:
      - The first row contains the backbone tokens.
      - The number of columns equals the number of backbone tokens.
      - The number of rows equals 1 (for the backbone) plus the maximum branch token length.
      - For each branch, its tokens are placed vertically in the column corresponding 
        to its attachment point (the last backbone "C" before the "("). 
      - Empty cells are filled with "_" (underscore).
    
    The output is written as plain text lines (tokens separated by a space) with no
    brackets or commas.
    """
    def process_segment(segment):
        """
        Process a segment (either backbone or branch substring) to:
          - Insert a single bond ("-") between consecutive "C"s if no explicit bond symbol is provided.
          - Preserve any explicit bond symbols: "-", "=" or "#".
        
        Returns a list of tokens.
        """
        tokens = []
        i = 0
        while i < len(segment):
            ch = segment[i]
            if ch == "C":
                # If the previous token is a "C", insert "-" before adding a new "C".
                if tokens and tokens[-1] == "C":
                    tokens.append("-")
                tokens.append("C")
                i += 1
            elif ch in ["-", "=", "#"]:
                tokens.append(ch)
                i += 1
            else:
                tokens.append(ch)
                i += 1
        return tokens

    # --- Step 1: Parse the SMILES string into backbone tokens and branches ---
    backbone_tokens = []  # tokens for the backbone (outside any parentheses)
    branches = []         # list of tuples: (attachment_column, branch_tokens)
    i = 0
    while i < len(smiles):
        if smiles[i] == '(':
            # A branch starts.
            # The branch attaches at the last "C" in backbone_tokens.
            attach_index = None
            for j in range(len(backbone_tokens) - 1, -1, -1):
                if backbone_tokens[j] == "C":
                    attach_index = j
                    break
            # Now, collect the branch substring until the matching ')'
            i += 1  # skip '('
            branch_str = ""
            paren_count = 1
            while i < len(smiles) and paren_count > 0:
                if smiles[i] == '(':
                    paren_count += 1
                elif smiles[i] == ')':
                    paren_count -= 1
                    if paren_count == 0:
                        i += 1  # skip the closing ')'
                        break
                if paren_count > 0:
                    branch_str += smiles[i]
                i += 1
            # Process the branch substring to get tokens.
            branch_tokens = process_segment(branch_str)
            # Insert a connector bond ("-") at the beginning if missing,
            # unless the branch explicitly starts with "=" or "#".
            if branch_tokens and branch_tokens[0] not in ["-", "=", "#"]:
                branch_tokens = ["-"] + branch_tokens
            # Save the branch tokens along with its attachment column.
            if attach_index is not None:
                branches.append((attach_index, branch_tokens))
        else:
            # Process characters outside parentheses (backbone).
            if smiles[i] == "C":
                if backbone_tokens and backbone_tokens[-1] == "C":
                    backbone_tokens.append("-")
                backbone_tokens.append("C")
            elif smiles[i] in ["-", "=", "#"]:
                backbone_tokens.append(smiles[i])
            else:
                backbone_tokens.append(smiles[i])
            i += 1

    # --- Step 2: Determine grid dimensions ---
    cols = len(backbone_tokens)
    max_branch_height = 0
    for attach_index, branch_tokens in branches:
        if len(branch_tokens) > max_branch_height:
            max_branch_height = len(branch_tokens)
    rows = 1 + max_branch_height  # 1 for the backbone row + branch height (0 if no branch)

    # --- Step 3: Build the grid ---
    # Initialize grid with "_" tokens.
    grid = [ ["_" for _ in range(cols)] for _ in range(rows) ]
    # Place the backbone in row 0.
    for c in range(cols):
        grid[0][c] = backbone_tokens[c]
    # For each branch, place its tokens vertically in the appropriate column.
    # If multiple branches attach at the same column, concatenate their tokens (separated by '/').
    branch_dict = {}  # maps column index -> list of branch token lists
    for attach_index, branch_tokens in branches:
        if attach_index not in branch_dict:
            branch_dict[attach_index] = []
        branch_dict[attach_index].append(branch_tokens)
    # For each attachment column, merge the branches row by row.
    for col, branch_token_lists in branch_dict.items():
        for r in range(max_branch_height):
            cell_tokens = []
            for tokens in branch_token_lists:
                if r < len(tokens):
                    cell_tokens.append(tokens[r])
            if cell_tokens:
                grid[r+1][col] = "/".join(cell_tokens)

    # --- Step 4: Write the grid to a text file ---
    with open(f"{compound_name}_atom.txt", "w") as f:
        for row in grid:
            line = " ".join(row)
            f.write(line + "\n")
            
def remove_special_characters(smiles: str) -> str:
    """
    Removes the characters 'H', '@', '[', ']', '+', and '-' from a given string.

    Parameters:
        smiles (str): The input string.

    Returns:
        str: The modified string with specified characters removed.
    """
    return ''.join(char for char in smiles if char not in {'H', '@', '[', ']', '+', '-'})


def transform_smiles(smiles: str) -> str:
    """
    Transforms a SMILES string by:
      - Replacing all capitalized letters that are not 'C' with 'C'.
      - Removing lowercase 'r' and 'l'.
      - Replacing all other lowercase letters with 'C'.

    Parameters:
        smiles (str): The input SMILES string.

    Returns:
        str: The transformed SMILES string.
    """
    transformed = []
    
    for char in smiles:
        if char.isupper():
            # Change all uppercase letters to 'C' if they are not already 'C'
            transformed.append('C' if char != 'C' else char)
        elif char.islower():
            # Remove 'r' and 'l', replace other lowercase with 'C'
            if char not in {'r', 'l'}:
                transformed.append('C')
        else:
            # Keep non-alphabetic characters unchanged
            transformed.append(char)
    
    return ''.join(transformed)