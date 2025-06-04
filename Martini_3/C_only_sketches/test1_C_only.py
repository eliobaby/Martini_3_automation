def get_c_only_dict():
    c_flipped = {}

    # ===================================================
    # 1. FULL BEADS (no continuation groups)
    # ===================================================

    # 2-carbon full beads:
    c_flipped["TF1"] = ["CC"]  # all single bonds
    c_flipped["TF2"] = ["C=C"]  # internal double → flag b

    # 3-carbon full beads:
    c_flipped["SF1"] = ["CCC"]
    c_flipped["SF2"] = ["C=CC"]   # internal double on left (canonical)
    c_flipped["SF3"] = ["C=C=C"]  # two internal doubles

    # 4-carbon full beads:
    c_flipped["LF1"] = ["CCCC"]
    c_flipped["LF2"] = ["C=CCC"]    # double at first bond → flag b
    c_flipped["LF3"] = ["CC=CC"]    # double in the middle → flag b
    c_flipped["LF4"] = ["C=CC=C"]   # two internal doubles → flag b
    c_flipped["LF5"] = ["C=C=CC"]   # alternate arrangement → flag b
    c_flipped["LF6"] = ["C=C=C=C"]  # all bonds double → flag b

    # Branched full beads (backbone not a straight line):
    # 4-carbon branched:
    c_flipped["LU1"] = ["CC(C)C"]     # linear 4-carbon branched variant (canonical: branch on second C)
    c_flipped["LU2"] = ["C=C(C)C"]    # with an internal double
    c_flipped["LU3"] = ["CC(=C)C"]

    # ===================================================
    # 2. EDGE BEADS (one continuation group at the right)
    # Pattern: main chain followed by a trailing continuation.
    # In the value string, if the trailing continuation is double ("++"), add flag d.
    # Also, if an internal double is present in the main chain, add b.
    # ===================================================

    # 2-carbon edge beads:
    c_flipped["TE1"]   = ["CC+"]         # single continuation
    c_flipped["TE2"]  = ["C=C+"]        # internal double → flag b
    c_flipped["TEd1"]  = ["CC++"]        # double continuation → flag d
    c_flipped["TEd2"] = ["C=C++"]       # both internal double and double continuation

    # 3-carbon edge beads:
    c_flipped["SE1"]    = ["CCC+"]
    c_flipped["SE2"]  = ["C=CC+"]      # internal double on left (C=CC+)
    c_flipped["SE3"]  = ["CC=C+"]      # internal double on right (CC=C+)
    c_flipped["SE4"]  = ["C=C=C+"]
    c_flipped["SEd1"]   = ["CCC++"]
    c_flipped["SEd2"] = ["C=CC++"]
    c_flipped["SEd3"] = ["CC=C++"]
    c_flipped["SEd4"] = ["C=C=C++"]

    # 4-carbon edge beads:
    c_flipped["LE1"]   = ["CCCC+"]
    c_flipped["LE2"]  = ["C=CCC+"]    # internal double → flag b
    c_flipped["LE3"]  = ["CC=CC+"]
    c_flipped["LE4"]  = ["CCC=C+"]
    c_flipped["LE5"]  = ["C=CC=C+"]   # two internal doubles → flag b
    c_flipped["LE6"]  = ["CC=C=C+"]
    c_flipped["LE7"]  = ["C=C=CC+"]
    c_flipped["LE8"]  = ["C=C=C=C+"]
    c_flipped["LEd1"] = ["CCCC++"]    # double continuation → flag d
    c_flipped["LEd2"] = ["C=CCC++"]
    c_flipped["LEd3"] = ["CC=CC++"]
    c_flipped["LEd4"] = ["CCC=C++"]
    c_flipped["LEd5"] = ["C=CC=C++"]
    c_flipped["LEd6"] = ["CC=C=C++"]
    c_flipped["LEd7"] = ["C=C=CC++"]
    c_flipped["LEd8"] = ["C=C=C=C++"]

    # ===================================================
    # 3. CONTINUATION BEADS (continuation groups at both ends)
    # Pattern: a leading continuation group, the main chain, then a trailing continuation.
    # For asymmetric cases (if one side is double), we use canonical form: left always "+" and right "++".
    # In the value, add flag d if the right continuation is double.
    # Internal double → add b.
    # ===================================================

    # 2-carbon continuation beads:
    c_flipped["TC1"]   = ["+CC+"]
    c_flipped["TCm1"]  = ["+CC++"]
    c_flipped["TCd1"]  = ["++CC++"]
    c_flipped["TC2"]  = ["+C=C+"]
    c_flipped["TCm2"] = ["+C=C++"]
    c_flipped["TCd2"] = ["++C=C++"]

    # 3-carbon continuation beads:
    c_flipped["SC1"]    = ["+CCC+"]
    c_flipped["SCm1"]   = ["+CCC++"]
    c_flipped["SCd1"]   = ["++CCC++"]
    c_flipped["SC2"]  = ["+C=CC+"]     # canonical: internal double on left (key: +C=CC+)
    c_flipped["SCm2"] = ["+C=CC++"]
    c_flipped["SCm3"] = ["+CC=C++"]
    c_flipped["SCd2"] = ["++C=CC++"]
    c_flipped["SC3"]  = ["+C=C=C+"]
    c_flipped["SCm4"] = ["+C=C=C++"]
    c_flipped["SCd3"] = ["++C=C=C++"]

    # 4-carbon continuation beads:
    c_flipped["LC1"]    = ["+CCCC+"]
    c_flipped["LCm1"]   = ["+CCCC++"]
    c_flipped["LCd1"]   = ["++CCCC++"]
    c_flipped["LC2"]  = ["+C=CCC+"]    # internal double → flag b
    c_flipped["LC3"]  = ["+CC=CC+"]
    c_flipped["LC4"]  = ["+C=CC=C+"]
    c_flipped["LC5"]  = ["+C=C=CC+"]
    c_flipped["LC6"]  = ["+C=C=C=C+"]
    c_flipped["LCm2"] = ["+C=CCC++"]
    c_flipped["LCm3"] = ["+CC=CC++"]
    c_flipped["LCm4"] = ["+CCC=C++"]
    c_flipped["LCm5"] = ["+C=CC=C++"]
    c_flipped["LCm6"] = ["+C=C=CC++"]
    c_flipped["LCm7"] = ["+CC=C=C++"]
    c_flipped["LCm8"] = ["+C=C=C=C++"]
    c_flipped["LCd2"] = ["++C=CCC++"]
    c_flipped["LCd3"] = ["++CC=CC++"]
    c_flipped["LCd4"] = ["++C=CC=C++"]
    c_flipped["LCd5"] = ["++C=C=CC++"]
    c_flipped["LCd6"] = ["++C=C=C=C++"]

    # ===================================================
    # 4. SIDECHAIN NODE BEADS (branch points with 3 connections)
    # For sidechain attachments we always use 2-atom beads, where one atom is in the backbone
    # and the other is the first atom of the sidechain.
    # Optional leading "+" indicates a left (backbone) continuation.
    # Optional trailing continuation ("+" for single, "++" for double) indicates right (backbone) continuation.
    # ===================================================

    # with all continuation (middle of the backbone)
    # CODE MID n-n
    c_flipped["TY1"]   = ["+C(C+)+"]      # All single bonds; left continuation implicit; right continuation single.
    c_flipped["TYm"]  = ["+C(C+)++"]     # Right continuation is double.
    c_flipped["TY2"]  = ["+C(=C+)+"]     # Bond to sidechain is double.
    c_flipped["TYs1"] = ["+C(C++)+"]     # All single bonds; left continuation implicit; right continuation single.
    c_flipped["TYms"] = ["+C(C++)++"]    # Right continuation is double.
    c_flipped["TYs2"] = ["+C(=C++)+"]    # Bond to sidechain is double.

    # ===================================================
    # Backbone with sidechain length 1 (continuation on both ends):
    # CODE MID 1-n
    c_flipped["TV1"]  = ["+C(C)+"]    # all single
    c_flipped["TV2"] = ["+C(=C)+"]   # internal double → flag b
    c_flipped["TVm"] = ["+C(C)++"]   # right continuation double → flag d

    # Backbone with sidechain length 2 (continuation on both ends):
    # CODE MID 2-n
    c_flipped["SV1"]   = ["+C(CC)+"]    # all single
    c_flipped["SV2"] = ["+C(C=C)+"]   # internal double → flag b
    c_flipped["SV3"] = ["+C(=CC)+"]
    c_flipped["SV4"] = ["+C(=C=C)+"]
    c_flipped["SVm1"]  = ["+C(CC)++"]
    c_flipped["SVm2"] = ["+C(C=C)++"]
    
    # Backbone with sidechain length 1 and an extra node on the backbone
    # CODE MID 1-n+1
    c_flipped["SW1"]  = ["+C(C)C+"]    # all single
    c_flipped["SWr"]  = ["++C(C)C+"]   
    c_flipped["SW2"] = ["+C(=C)C+"]   # internal double → flag b
    c_flipped["SW3"] = ["+C(C)=C+"]   # right continuation double → flag d
    c_flipped["SWm1"]  = ["+C(C)C++"]    # all single
    c_flipped["SWd"]  = ["++C(C)C++"]   
    c_flipped["SWm2"] = ["+C(=C)C++"]   # internal double → flag b
    c_flipped["SWm3"] = ["+C(C)=C++"]   # right continuation double → flag d
    
    # Backbone with sidechain length 2 and an extra node on the backbone
    # CODE MID n-n+1
    c_flipped["TW1"]  = ["+C(+)C+"]    # all single
    c_flipped["TWr"]  = ["++C(+)C+"]   
    c_flipped["TWs"] = ["+C(++)C+"]   # internal double → flag b
    c_flipped["TW2"] = ["+C(+)=C+"]   # right continuation double → flag d
    c_flipped["TWm1"]  = ["+C(+)C++"]    # all single
    c_flipped["TWd"]  = ["++C(+)C++"]   
    c_flipped["TWms"] = ["+C(++)C++"]   # internal double → flag b
    c_flipped["TWm2"] = ["+C(+)=C++"]   # right continuation double → flag d

####################################################### FINISH
    # ===================================================
    # 5. TRIPLE BOND BEADS (using '#' notation; unchanged)
    # ===================================================
    c_flipped["TX1"] = ["C#C"] 
    c_flipped["TX2"] = ["C#C+"] 
    c_flipped["SX4"] = ["C#CC+"] 
    c_flipped["SX5"] = ["C#CC++"] 
    c_flipped["SX1"] = ["C#CC"] 
    c_flipped["TX3"] = ["+C#C+"] 
    c_flipped["SX2"] = ["+C#CC+"] #why is this here (very important)
    c_flipped["SXm"] = ["+C#CC++"] #why is this here (very important)
    c_flipped["SX3"] = ["CC#C+"] 
    c_flipped["LX2"] = ["CC#CC+"]
    c_flipped["LX3"] = ["CC#CC++"] 
    c_flipped["LX1"] = ["CC#CC"] 

    # ===================================================
    # 6. EXTREME EDGE CASES
    # ===================================================
    # Sidechain with length 1, 1 away from the end of a backbone
    # CODE LEFT 1-1 CODE RIGHT 1-1 
    c_flipped["SZ1"]  = ["CC(C)+"]
    c_flipped["SZ2"] = ["C=C(C)+"]
    c_flipped["SZ3"] = ["CC(=C)+"]
    c_flipped["SZd"] = ["CC(C)++"]
    # The next node is lonely:
    # CODE LEFT 1-1+1
    c_flipped["LZ1"]  = ["CC(C)C+"]
    c_flipped["LZ2"] = ["C=C(C)C+"]
    c_flipped["LZ3"] = ["CC(=C)C+"]
    c_flipped["LZ4"] = ["CC(C)=C+"]
    c_flipped["LZd1"]  = ["CC(C)C++"]
    c_flipped["LZd2"] = ["C=C(C)C++"]
    c_flipped["LZd3"] = ["CC(=C)C++"]
    c_flipped["LZd4"] = ["CC(C)=C++"]
    
    # Sidechain with length n, 1 away from the end of a backbone
    # CODE LEFT 1-n CODE RIGHT n-1
    c_flipped["TH1"] = ["CC(+)+"]
    c_flipped["TH2"] = ["C=C(+)+"]
    c_flipped["THs"] = ["CC(++)+"]
    c_flipped["THd"] = ["CC(+)++"]
    # The next node is lonely:
    # CODE LEFT 1-n+1
    c_flipped["SH1"] = ["CC(+)C+"]
    c_flipped["SH2"] = ["C=C(+)C+"]
    c_flipped["SHs"] = ["CC(++)C+"]
    c_flipped["SH3"] = ["CC(+)=C+"]
    c_flipped["SHd1"] = ["CC(+)C++"]
    c_flipped["SHd2"] = ["C=C(+)C++"]
    c_flipped["SHds"] = ["CC(++)C++"]
    c_flipped["SHd3"] = ["CC(+)=C++"]
    # note to self: n could be 0 for pesky triple bonds
    c_flipped["OF"] = ["C"]

    return c_flipped






