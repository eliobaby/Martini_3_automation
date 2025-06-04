import random
import string

def generate_random_string(length=6):
    """
    Generates a random string of a given length using all printable characters.
    
    Parameters:
        length (int): The length of the generated string (default is 6).
    
    Returns:
        str: A random string of the specified length.
    """
    allowed_chars = string.printable.strip()  # Excludes whitespace characters
    return ''.join(random.choices(allowed_chars, k=length))

def fix_matrix_and_bead(matrix_full, bead_empty, bead, starting_location):
    """
    Updates matrix_full and bead_empty based on the provided bead string and starting location.
    
    Input:
      - matrix_full (3d list)
      - bead_empty (3d list)
      - bead (string)
      - starting_location (tuple (a, b))
      
    The function uses a random string (generated via generate_random_string) and applies different
    modifications to the matrices based on the following cases:
    
      Case 1: If 'OF' is found in bead:
          - Set bead_empty[a][b][0] = bead and bead_empty[a][b][1] = random_string.
          - Set matrix_full[a][b][0] = 'X'.
          
      Case 2: If any of 'TF', 'TE', 'TC', 'TW', 'TX', or 'TH' is found in bead:
          - If a != 0:
                * Update bead_empty[a][b] and bead_empty[a+1][b] with bead and random_string.
                * Set matrix_full[a][b][0] and matrix_full[a+1][b][0] to 'X'.
          - Else:
                * Update bead_empty[a][b] and bead_empty[a][b+1] similarly.
                * Set matrix_full[a][b][0] and matrix_full[a][b+1][0] to 'X'.
          
      Case 3: If any of 'SF', 'SE', 'SC', 'SH', or 'SX' is found in bead:
          - If a != 0:
                * Update bead_empty at rows a, a+1, a+2 (same column b) with bead and random_string.
                * Set matrix_full at these positions in column b to 'X'.
          - Else:
                * Update bead_empty at column indices b, b+1, b+2 (same row a) with bead and random_string.
                * Set matrix_full at these positions to 'X'.
          
      Case 4: If any of 'LF', 'LE', 'LC', or 'LX' is found in bead:
          - If a != 0:
                * Update bead_empty at rows a, a+1, a+2, a+3 (column b) with bead and random_string.
                * Set matrix_full at these positions to 'X'.
          - Else:
                * Update bead_empty at columns b, b+1, b+2, b+3 (row a) with bead and random_string.
                * Set matrix_full at these positions to 'X'.
          
      Case 5: If either 'LU' or 'LZ' is found in bead:
                * Update bead_empty at [a][b], [a][b+1], [a][b+2] and at [a+1][b+1] with bead and random_string.
                * Set the corresponding positions in matrix_full to 'X'.
          
      Case 6: If either 'TY' or 'TV' is found in bead:
                * Update bead_empty at [a][b] and [a+1][b] with bead and random_string.
                * Set matrix_full at these positions to 'X'.
          
      Case 7: If 'SW' is found in bead:
                * Update bead_empty at [a][b], [a+1][b], and [a][b+1] with bead and random_string.
                * Set matrix_full at these positions to 'X'.
          
      Case 8: If 'SV' is found in bead:
                * Update bead_empty at [a][b], [a+1][b], and [a+2][b] with bead and random_string.
                * Set matrix_full at these positions to 'X'.
          
      Case 9: If 'SZ' is found in bead:
                * Update bead_empty at [a][b], [a+1][b], and [a+1][b+1] with bead and random_string.
                * Set matrix_full at these positions to 'X'.
      
    Returns:
      tuple: (matrix_full, bead_empty) after modifications.
    """
    a, b = starting_location
    rand_str = generate_random_string()

    # Case 1: 'OF'
    if 'OF' in bead:
        bead_empty[a][b][0] = bead
        bead_empty[a][b][1] = rand_str
        matrix_full[a][b][0] = 'X'

    # Case 2: 'TF', 'TE', 'TC', 'TW', 'TX', 'TH'
    elif any(sub in bead for sub in ['TF', 'TE', 'TC', 'TW', 'TX', 'TH']):
        if a != 0:
            bead_empty[a][b][0] = bead
            bead_empty[a][b][1] = rand_str
            bead_empty[a+1][b][0] = bead
            bead_empty[a+1][b][1] = rand_str
            matrix_full[a][b][0] = 'X'
            matrix_full[a+1][b][0] = 'X'
        else:
            bead_empty[a][b][0] = bead
            bead_empty[a][b][1] = rand_str
            bead_empty[a][b+1][0] = bead
            bead_empty[a][b+1][1] = rand_str
            matrix_full[a][b][0] = 'X'
            matrix_full[a][b+1][0] = 'X'

    # Case 3: 'SF', 'SE', 'SC', 'SH', 'SX'
    elif any(sub in bead for sub in ['SF', 'SE', 'SC', 'SH', 'SX']):
        if a != 0:
            bead_empty[a][b][0]     = bead; bead_empty[a][b][1]     = rand_str
            bead_empty[a+1][b][0]   = bead; bead_empty[a+1][b][1]   = rand_str
            bead_empty[a+2][b][0]   = bead; bead_empty[a+2][b][1]   = rand_str
            matrix_full[a][b][0]    = 'X'
            matrix_full[a+1][b][0]  = 'X'
            matrix_full[a+2][b][0]  = 'X'
        else:
            bead_empty[a][b][0]     = bead; bead_empty[a][b][1]     = rand_str
            bead_empty[a][b+1][0]   = bead; bead_empty[a][b+1][1]   = rand_str
            bead_empty[a][b+2][0]   = bead; bead_empty[a][b+2][1]   = rand_str
            matrix_full[a][b][0]    = 'X'
            matrix_full[a][b+1][0]  = 'X'
            matrix_full[a][b+2][0]  = 'X'

    # Case 4: 'LF', 'LE', 'LC', 'LX'
    elif any(sub in bead for sub in ['LF', 'LE', 'LC', 'LX']):
        if a != 0:
            bead_empty[a][b][0]     = bead; bead_empty[a][b][1]     = rand_str
            bead_empty[a+1][b][0]   = bead; bead_empty[a+1][b][1]   = rand_str
            bead_empty[a+2][b][0]   = bead; bead_empty[a+2][b][1]   = rand_str
            bead_empty[a+3][b][0]   = bead; bead_empty[a+3][b][1]   = rand_str
            matrix_full[a][b][0]    = 'X'
            matrix_full[a+1][b][0]  = 'X'
            matrix_full[a+2][b][0]  = 'X'
            matrix_full[a+3][b][0]  = 'X'
        else:
            bead_empty[a][b][0]     = bead; bead_empty[a][b][1]     = rand_str
            bead_empty[a][b+1][0]   = bead; bead_empty[a][b+1][1]   = rand_str
            bead_empty[a][b+2][0]   = bead; bead_empty[a][b+2][1]   = rand_str
            bead_empty[a][b+3][0]   = bead; bead_empty[a][b+3][1]   = rand_str
            matrix_full[a][b][0]    = 'X'
            matrix_full[a][b+1][0]  = 'X'
            matrix_full[a][b+2][0]  = 'X'
            matrix_full[a][b+3][0]  = 'X'

    # Case 5: 'LU' or 'LZ'
    elif any(sub in bead for sub in ['LU', 'LZ']):
        bead_empty[a][b][0]      = bead; bead_empty[a][b][1]      = rand_str
        bead_empty[a][b+1][0]    = bead; bead_empty[a][b+1][1]    = rand_str
        bead_empty[a][b+2][0]    = bead; bead_empty[a][b+2][1]    = rand_str
        bead_empty[a+1][b+1][0]  = bead; bead_empty[a+1][b+1][1]  = rand_str
        matrix_full[a][b][0]     = 'X'
        matrix_full[a][b+1][0]   = 'X'
        matrix_full[a][b+2][0]   = 'X'
        matrix_full[a+1][b+1][0] = 'X'

    # Case 6: 'TY' or 'TV'
    elif any(sub in bead for sub in ['TY', 'TV']):
        bead_empty[a][b][0] = bead; bead_empty[a][b][1] = rand_str
        bead_empty[a+1][b][0] = bead; bead_empty[a+1][b][1] = rand_str
        matrix_full[a][b][0] = 'X'
        matrix_full[a+1][b][0] = 'X'

    # Case 7: 'SW'
    elif "SW" in bead:
        bead_empty[a][b][0]   = bead; bead_empty[a][b][1]   = rand_str
        bead_empty[a+1][b][0] = bead; bead_empty[a+1][b][1] = rand_str
        bead_empty[a][b+1][0] = bead; bead_empty[a][b+1][1] = rand_str
        matrix_full[a][b][0]  = 'X'
        matrix_full[a+1][b][0]= 'X'
        matrix_full[a][b+1][0]= 'X'

    # Case 8: 'SV'
    elif "SV" in bead:
        bead_empty[a][b][0]    = bead; bead_empty[a][b][1]    = rand_str
        bead_empty[a+1][b][0]  = bead; bead_empty[a+1][b][1]  = rand_str
        bead_empty[a+2][b][0]  = bead; bead_empty[a+2][b][1]  = rand_str
        matrix_full[a][b][0]   = 'X'
        matrix_full[a+1][b][0] = 'X'
        matrix_full[a+2][b][0] = 'X'

    # Case 9: 'SZ'
    elif "SZ" in bead:
        bead_empty[a][b][0]     = bead; bead_empty[a][b][1]     = rand_str
        bead_empty[a+1][b][0]   = bead; bead_empty[a+1][b][1]   = rand_str
        bead_empty[a+1][b+1][0] = bead; bead_empty[a+1][b+1][1] = rand_str
        matrix_full[a][b][0]    = 'X'
        matrix_full[a+1][b][0]  = 'X'
        matrix_full[a+1][b+1][0]= 'X'

    # If no case is matched, do nothing (or you could handle it differently)
    return matrix_full, bead_empty
    
def full_bead_assignment(matrix_full, bead_empty):
    """
    Processes the full bead section (2-4 nodes) by checking the dimensions of matrix_full
    and inspecting specific cell values, then calls fix_matrix_and_bead with the appropriate
    bead string and starting location.

    Parameters:
      matrix_full (list): A 3D list representing the full matrix.
      bead_empty (list): A 3D list of the same dimensions as matrix_full.

    Returns:
      tuple: (matrix_full, bead_empty) after modification.
    """
    num_rows = len(matrix_full)
    num_cols = len(matrix_full[0]) if num_rows > 0 else 0

    # Case 1: Size is 1x2
    if num_rows == 1 and num_cols == 2:
        if matrix_full[0][0][3] == 1:
            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TF1', (0, 0))
        else:
            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TF2', (0, 0))
    
    # Case 2: Size is 1x3
    elif num_rows == 1 and num_cols == 3:
        if matrix_full[0][0][3] == 1 and matrix_full[0][1][3] == 1:
            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SF1', (0, 0))
        elif matrix_full[0][0][3] == 1 or matrix_full[0][1][3] == 1:
            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SF2', (0, 0))
        else:
            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SF3', (0, 0))
    
    # Case 3: Size is 1x4
    elif num_rows == 1 and num_cols == 4:
        a = matrix_full[0][0][3]
        b = matrix_full[0][1][3]
        c = matrix_full[0][2][3]
        if a == 1 and b == 1 and c == 1:
            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'LF1', (0, 0))
        elif (a == 1 and b == 1) or (b == 1 and c == 1):
            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'LF2', (0, 0))
        elif a == 1 and c == 1:
            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'LF3', (0, 0))
        elif a == 1 or c == 1:
            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'LF5', (0, 0))
        elif b == 1:
            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'LF4', (0, 0))
        else:
            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'LF6', (0, 0))
    
    # Case 4: Size is 2x3
    elif num_rows == 2 and num_cols == 3:
        # Here we use the second row (index 1) of the first column's bead info.
        a = matrix_full[0][1][1]
        b = matrix_full[0][1][2]
        c = matrix_full[0][1][3]
        if a == 2 or c == 2:
            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'LU2', (0, 0))
        elif b == 2:
            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'LU3', (0, 0))
        else:
            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'LU1', (0, 0))
    
    return matrix_full, bead_empty

def triple_bond_assignment(matrix_full, bead_empty):
    """
    Assigns triple bond beads to the given matrices based on conditions in matrix_full.
    
    Input:
      - matrix_full: a 3D list (for example, backbone and sidechain info)
      - bead_empty: a 3D list of the same dimensions as matrix_full
    Returns:
      (matrix_full, bead_empty) after modifications.
      
    This function uses the following helper function:
      fix_matrix_and_bead(matrix_full, bead_empty, bead, starting_location)
    which is assumed to update both matrices based on the bead string and starting location.
    """
    num_cols = len(matrix_full[0])
    
    # Process row index 0 (backbone) first.
    for col in range(num_cols):
        cell = matrix_full[0][col]
        # If cell has a 'C' in index 0 and a 3 in index 3:
        if cell[0] == 'C' and cell[3] == 3:
            # Case: 2D size is 4 (i.e. 1x4 backbone)
            if num_cols == 4:
                if col == 1:
                    matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'LX1', (0, 0))                    
                elif col == 0:
                    matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX2', (0, 0))                    
                elif col == 2:
                    matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX2', (0, 2))
                    
            # Case: 2D size is 2 (1x2)
            elif num_cols == 2:
                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX1', (0, 0))                
            # Case: 2D size is 3 (1x3)
            elif num_cols == 3:
                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX1', (0, 0))                
            else:
                # Else if the number 3 is found in index 3 for other column positions:
                if col == 0:
                    if len(matrix_full) > 1 and len(matrix_full[1]) > 3:
                        if matrix_full[1][2][0] == 'X' and matrix_full[1][3][0] == 'C':
                            if len(matrix_full[0]) > 2 and matrix_full[0][2][3] == 1:
                                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX4', (0, 0))                                
                            elif len(matrix_full[0]) > 2 and matrix_full[0][2][3] == 2:
                                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX5', (0, 0))                                
                            else:
                                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX2', (0, 0))
                        else:
                            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX2', (0, 0))    
                    else:
                        matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX2', (0, 0))    
                elif col == num_cols - 2:  # end column - 1
                    matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX2', (0, col))
                elif col == 1:
                    if len(matrix_full) > 1 and len(matrix_full[1]) > 4:
                        if matrix_full[1][3][0] == 'X' and matrix_full[1][4][0] == 'C':
                            if matrix_full[0][3][3] == 1:
                                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'LX2', (0, 0))
                            elif matrix_full[0][3][3] == 2:
                                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'LX3', (0, 0))
                            else:
                                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX3', (0, 0))
                        else:
                            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX3', (0, 0))
                    else:
                        matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX3', (0, 0))
                elif col == num_cols - 3:  # end column - 2
                    matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX3', (0, col))
                else:
                    # Otherwise, use current column as reference:
                    col_index = col
                    if len(matrix_full[1]) > col_index + 3:
                        if (matrix_full[1][col_index+2][0] == 'X' and 
                            matrix_full[1][col_index+3][0] == 'C'):
                            if matrix_full[0][col_index+2][3] == 1:
                                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX2', (0, col_index))
                            elif matrix_full[0][col_index+2][3] == 2:
                                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SXm', (0, col_index))
                            else:
                                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX3', (0, col_index))
                        else:
                            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX3', (0, col_index))
                    else:
                        matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX3', (0, col_index))      

    # Process row index 1 (sidechain region)
    if len(matrix_full) > 1:
        for col in range(len(matrix_full[1])):
            if matrix_full[1][col][0] == 'X':
                continue
            else:
                for row in range(1, len(matrix_full)):
                    if matrix_full[row][col][0] == 'C' and matrix_full[row][col][3] == 3:
                        # First branch: check if row+2 and row+3 exist and have 'C'
                        if row+2 < len(matrix_full):
                            if row+3 < len(matrix_full):
                                if (matrix_full[row+2][col][0] == 'C' and 
                                    matrix_full[row+3][col][0] == 'C'):
                                    if row == 2:
                                        valid_condition = False
                                        if (row-1 >= 0 and col+1 < len(matrix_full[row-1]) and matrix_full[row-1][col+1][0] == 'X' and 
                                            col+2 < len(matrix_full[row-1]) and matrix_full[row-1][col+2][0] == 'C'):
                                            valid_condition = True
                                        elif (row-2 >= 0 and (col+2 >= len(matrix_full[row-2]) or matrix_full[row-2][col+2][0] is None)):
                                            valid_condition = True
                                        if valid_condition:
                                            if matrix_full[row-2][col][3] == 1:
                                                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX2', (row-1, col))
                                            elif matrix_full[row-2][col][3] == 2:
                                                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SXm', (row-1, col))
                                        else:
                                            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX3', (row, col))
                                    else:
                                        matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX3', (row, col))
                                elif matrix_full[row+2][col][0] == 'C':
                                    # Second branch:
                                    if row == 2:
                                        valid_condition = False
                                        if (row-1 >= 0 and col+1 < len(matrix_full[row-1]) and matrix_full[row-1][col+1][0] == 'X' and 
                                            col+2 < len(matrix_full[row-1]) and matrix_full[row-1][col+2][0] == 'C'):
                                            valid_condition = True
                                        elif (row-2 >= 0 and (col+2 >= len(matrix_full[row-2]) or matrix_full[row-2][col+2][0] is None)):
                                            valid_condition = True
                                        if valid_condition:
                                            if matrix_full[row-2][col][3] == 1:
                                                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'LX2', (row-1, col))
                                            elif matrix_full[row-2][col][3] == 2:
                                                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'LX3', (row-1, col))
                                        else:
                                            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX3', (row, col))
                                    else:
                                        matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX3', (row, col))
                                # Third branch:
                                else:
                                    if row == 2:
                                        valid_condition = False
                                        if (row-1 >= 0 and col+1 < len(matrix_full[row-1]) and matrix_full[row-1][col+1][0] == 'X' and 
                                            col+2 < len(matrix_full[row-1]) and matrix_full[row-1][col+2][0] == 'C'):
                                            valid_condition = True
                                        elif (row-2 >= 0 and (col+2 >= len(matrix_full[row-2]) or matrix_full[row-2][col+2][0] is None)):
                                            valid_condition = True
                                        if valid_condition:
                                            if matrix_full[row-2][col][3] == 1:
                                                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX4', (row-1, col))
                                            elif matrix_full[row-2][col][3] == 2:
                                                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX5', (row-1, col))
                                        else:
                                            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX2', (row, col))
                                    else:
                                        matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX2', (row, col))
                            elif matrix_full[row+2][col][0] == 'C':
                                # Second branch:
                                if row == 2:
                                    valid_condition = False
                                    if (row-1 >= 0 and col+1 < len(matrix_full[row-1]) and matrix_full[row-1][col+1][0] == 'X' and 
                                        col+2 < len(matrix_full[row-1]) and matrix_full[row-1][col+2][0] == 'C'):
                                        valid_condition = True
                                    elif (row-2 >= 0 and (col+2 >= len(matrix_full[row-2]) or matrix_full[row-2][col+2][0] is None)):
                                        valid_condition = True
                                    if valid_condition:
                                        if matrix_full[row-2][col][3] == 1:
                                            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'LX2', (row-1, col))
                                        elif matrix_full[row-2][col][3] == 2:
                                            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'LX3', (row-1, col))
                                    else:
                                        matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX3', (row, col))
                                else:
                                    matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX3', (row, col))
                            # Third branch:
                            else:
                                if row == 2:
                                    valid_condition = False
                                    if (row-1 >= 0 and col+1 < len(matrix_full[row-1]) and matrix_full[row-1][col+1][0] == 'X' and 
                                        col+2 < len(matrix_full[row-1]) and matrix_full[row-1][col+2][0] == 'C'):
                                        valid_condition = True
                                    elif (row-2 >= 0 and (col+2 >= len(matrix_full[row-2]) or matrix_full[row-2][col+2][0] is None)):
                                        valid_condition = True
                                    if valid_condition:
                                        if matrix_full[row-2][col][3] == 1:
                                            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX4', (row-1, col))
                                        elif matrix_full[row-2][col][3] == 2:
                                            matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX5', (row-1, col))
                                    else:
                                        matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX2', (row, col))
                                else:
                                    matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX2', (row, col))
                        # Third branch:
                        else:
                            if row == 2:
                                valid_condition = False
                                if (row-1 >= 0 and col+1 < len(matrix_full[row-1]) and matrix_full[row-1][col+1][0] == 'X' and 
                                    col+2 < len(matrix_full[row-1]) and matrix_full[row-1][col+2][0] == 'C'):
                                    valid_condition = True
                                elif (row-2 >= 0 and (col+2 >= len(matrix_full[row-2]) or matrix_full[row-2][col+2][0] is None)):
                                    valid_condition = True
                                if valid_condition:
                                    if matrix_full[row-2][col][3] == 1:
                                        matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX4', (row-1, col))
                                    elif matrix_full[row-2][col][3] == 2:
                                        matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'SX5', (row-1, col))
                                else:
                                    matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX2', (row, col))
                            else:
                                matrix_full, bead_empty = fix_matrix_and_bead(matrix_full, bead_empty, 'TX2', (row, col))
    return matrix_full, bead_empty

