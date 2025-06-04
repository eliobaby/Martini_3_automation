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