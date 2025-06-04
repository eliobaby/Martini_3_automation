def comb_check(s: str) -> None:
    """
    Checks if the string s is a valid comb polymer.
    
    Allowed characters: 'C', '-', '=', '#', '(', ')'
    
    The function:
      - Increments a counter by 1 for '(' and decrements by 1 for ')'.
      - Immediately throws an error if the counter becomes any value other than 0 or 1.
      - Throws an error if:
            • s contains any character not in the allowed set,
            • s contains the substring ")(",
            • or after processing, the counter is not 0 or 1.
            
    If all checks pass, it prints:
      "This is a valid comb polymer and continues"
      
    Note: This function does not return a value.
    """
    
    allowed_chars = {'C', '-', '=', '#', '(', ')'}
    
    # Check for any invalid characters.
    for char in s:
        if char not in allowed_chars:
            raise ValueError("Polymer does not match the requirements.")
    
    # Check for the forbidden substring.
    if ")(" in s:
        raise ValueError("Polymer does not match the requirements.")
    
    # Trace the string and update count.
    count = 0
    for char in s:
        if char == '(':
            count += 1
        elif char == ')':
            count -= 1
        
        # Immediately error out if count is not 0 or 1.
        if count not in (0, 1):
            raise ValueError("Polymer does not match the requirements.")

    # If all checks pass, display the valid message.
    print("This is a valid comb polymer")
    
def matrix_size_calc(smiles: str) -> list:
    """
    Calculates the size of a matrix based on the given SMILES string and returns a 2D list.
    
    Process for determining the number of rows (row):
      - Initialize row to 1.
      - Traverse the string and start a temporary counter (row_temp) at 1.
      - When a '(' is encountered, begin counting: for every 'C' found while inside the branch, increment row_temp.
      - When a ')' is encountered, compare row_temp to the current row value; if row_temp is larger, update row.
        Then reset row_temp to 1.
    
    Process for determining the number of columns (column):
      - Initialize column to 0.
      - Traverse the string and, for every 'C' found that is not inside a branch (i.e. not between '(' and ')'),
        increment column.
    
    Finally, returns a 2D list of dimensions row x column where each element is a list: ['C', 0, 0, 0].
    """
    # Calculate the number of rows
    row = 1
    row_temp = 1
    inside_branch = False
    
    for char in smiles:
        if char == '(':
            inside_branch = True
            row_temp = 1  # reset temporary row count when branch starts
        elif char == ')':
            inside_branch = False
            if row_temp > row:
                row = row_temp
            row_temp = 1  # reset for next branch
        elif inside_branch:
            if char == 'C':
                row_temp += 1

    # Calculate the number of columns (only count 'C's outside branches)
    column = 0
    inside_branch = False
    for char in smiles:
        if char == '(':
            inside_branch = True
        elif char == ')':
            inside_branch = False
        elif char == 'C' and not inside_branch:
            column += 1

    # Create and return the matrix of size (row x column)
    # Each element of the matrix is a list: ['C', 0, 0, 0]
    matrix = [[['X', 0, 0, 0] for _ in range(column)] for _ in range(row)]
    return matrix

def backbone_maker(smiles: str, matrix_empty: list) -> list:
    """
    Processes a SMILES string and updates the provided matrix_empty to form the backbone of a comb polymer.
    
    Steps:
      1. Remove all parentheses and the characters enclosed by each matching pair. This extracts the main chain.
      2. Initialize a column counter at 0.
      3. Iterate over the cleaned (backbone) string:
         - For each index where the character is 'C':
             • Set matrix_empty[0][column][0] = 'C'
             • If a previous character exists:
                   - If it is 'C', set matrix_empty[0][column][1] = 1
                   - If it is '=', set matrix_empty[0][column][1] = 2
                   - If it is '#', set matrix_empty[0][column][1] = 3
             • If a next character exists:
                   - If it is 'C', set matrix_empty[0][column][3] = 1
                   - If it is '=', set matrix_empty[0][column][3] = 2
                   - If it is '#', set matrix_empty[0][column][3] = 3
             • Increment the column counter.
      4. Return the updated matrix_empty.
      
    Parameters:
      smiles (str): The input SMILES string.
      matrix_empty (list): A 2D list (with at least one row) pre-initialized to the correct size, where each element is a list in the form ['C', 0, 0, 0].
      
    Returns:
      list: The modified matrix_empty with backbone information inserted.
    """
    # Remove parentheses and everything inside them to extract the backbone.
    backbone = ""
    skip = 0  # Counter to track if we are inside parentheses.
    for char in smiles:
        if char == '(':
            skip += 1
        elif char == ')':
            if skip > 0:
                skip -= 1
        else:
            if skip == 0:
                backbone += char

    # Start column count at 0.
    column = 0
    # Traverse the backbone string.
    for i, char in enumerate(backbone):
        if char == 'C':
            # Set the current cell's first element to 'C'
            matrix_empty[0][column][0] = 'C'
            
            # Check previous character if it exists.
            if i > 0:
                prev_char = backbone[i - 1]
                if prev_char == 'C':
                    matrix_empty[0][column][1] = 1
                elif prev_char == '=':
                    matrix_empty[0][column][1] = 2
                elif prev_char == '#':
                    matrix_empty[0][column][1] = 3

            # Check next character if it exists.
            if i < len(backbone) - 1:
                next_char = backbone[i + 1]
                if next_char == 'C':
                    matrix_empty[0][column][3] = 1
                elif next_char == '=':
                    matrix_empty[0][column][3] = 2
                elif next_char == '#':
                    matrix_empty[0][column][3] = 3
            
            # Move to the next column.
            column += 1

    return matrix_empty

def sidechain_maker(smiles: str, matrix_backbone: list) -> list:
    """
    Processes a SMILES string and updates the provided matrix_backbone with sidechain information.
    
    The function works as follows:
      - It uses a column counter (col_count) starting at -1 and a row counter (row_count) starting at 0.
      - It iterates through the SMILES string while maintaining a flag (inside_branch) that is set True
        when a '(' is encountered and set False when a ')' is encountered.
      - For every character in the SMILES string:
          • If the character is 'C' and we are not inside a branch (flag is False), increment col_count.
          • If the character is '(':
                - Set the flag (inside_branch) to True.
                - Look at the next character (if it exists); if it is 'C', update matrix_backbone[0][col_count][2] to 1;
                  if it is '=', update that cell to 2.
          • If the character is 'C' and we are inside a branch (flag is True):
                - Increment row_count.
                - Update the current sidechain cell: matrix_backbone[row_count][col_count][0] is set to 'C'.
                - Check the previous character (if it exists):
                      - If it is 'C' or '(' then set matrix_backbone[row_count][col_count][1] = 1.
                      - If it is '=' then set it to 2.
                      - If it is '#' then set it to 3.
                - Check the next character (if it exists):
                      - If it is 'C' then set matrix_backbone[row_count][col_count][3] = 1.
                      - If it is '=' then set it to 2.
                      - If it is '#' then set it to 3.
          • If the character is ')':
                - Turn the flag off and reset row_count to 0.
                
      - Finally, the function returns the modified matrix_backbone.
      
    Parameters:
      smiles (str): The SMILES string containing the comb polymer notation.
      matrix_backbone (list): A pre-constructed 2D list (matrix) with dimensions sufficient to hold the
                              backbone (row 0) and sidechains (rows > 0) for each backbone carbon.
                              Each element is assumed to be a list of the form ['C', 0, 0, 0].
    
    Returns:
      list: The updated matrix_backbone with sidechain information filled in.
    """
    col_count = -1
    row_count = 0
    inside_branch = False
    n = len(smiles)
    
    for i, char in enumerate(smiles):
        if char == '(':
            # Start of sidechain: set flag and update backbone bond info for sidechain attachment.
            inside_branch = True
            # Update the backbone cell (row 0, current col) with bond type from next character, if available.
            if i + 1 < n:
                next_char = smiles[i + 1]
                if next_char == 'C':
                    matrix_backbone[0][col_count][2] = 1
                elif next_char == '=':
                    matrix_backbone[0][col_count][2] = 2
            continue  # move to next character

        if char == ')':
            # End of sidechain: reset flag and row counter.
            inside_branch = False
            row_count = 0
            continue
        
        if not inside_branch:
            # When not inside a branch, update backbone column count on encountering a 'C'.
            if char == 'C':
                col_count += 1
        else:
            # When inside a branch (sidechain region)
            if char == 'C':
                row_count += 1  # Move to the next row in the sidechain for this branch.
                # Set the cell's first element to 'C'
                matrix_backbone[row_count][col_count][0] = 'C'
                
                # Check the previous character in the SMILES string.
                if i - 1 >= 0:
                    prev_char = smiles[i - 1]
                    if prev_char == 'C' or prev_char == '(':
                        matrix_backbone[row_count][col_count][1] = 1
                    elif prev_char == '=':
                        matrix_backbone[row_count][col_count][1] = 2
                    elif prev_char == '#':
                        matrix_backbone[row_count][col_count][1] = 3
                
                # Check the next character in the SMILES string.
                if i + 1 < n:
                    next_char = smiles[i + 1]
                    if next_char == 'C':
                        matrix_backbone[row_count][col_count][3] = 1
                    elif next_char == '=':
                        matrix_backbone[row_count][col_count][3] = 2
                    elif next_char == '#':
                        matrix_backbone[row_count][col_count][3] = 3
    
    return matrix_backbone

def bead_empty_matrix(matrix_empty):
    """
    Creates a matrix with the same dimensions (rows and columns) as matrix_empty,
    where each element is an empty string.

    Parameters:
        matrix_empty (list): A 2D list representing the original matrix.

    Returns:
        list: A new 2D list with the same shape, filled with empty strings.
    """
    rows = len(matrix_empty)
    cols = len(matrix_empty[0])
    return [[['',''] for _ in range(cols)] for _ in range(rows)]