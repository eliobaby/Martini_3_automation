import random
import string
from typing import List, Dict, Any, Tuple, Optional
import warnings

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

def pick_bead_key(martini_dict: Dict[str, List[Any]],
                  element: str,
                  bond_order: Optional[int] = None,
                  kind: str = 'S') -> Optional[str]:
    """
    Pick a bead‐type key from martini_dict whose
      • key starts with `kind` ('S' or 'T'),
      • val[0] == 2,
      • and whose val[2] (the "(...)" part) matches our element + optional '='.

    element: e.g. 'O', 'N', 'S', 'CL', 'I', 'C', 'BR', F
    bond_order: if element == 'O' and bond_order == 2, we look for '(=O)' else '(O)'
    kind: 'S' for SN6*, 'SC6'…, 'SN6a', etc.; 'T' for TN6*, 'TC4'…, 'TX2', etc.
    """
    # build the pattern to look for inside the parentheses
    el = element.capitalize()  # so 'CL' -> 'Cl'
    if element == 'O' and bond_order == 2:
        pat = '(=O)'
    else:
        pat = f'({el})'

    for key, val in martini_dict.items():
        if not key.startswith(kind):
            continue
        if val[0] != 2:
            continue
        # val[2] is something like 'CC(O)' or 'C(Cl)' or 'CC(=O)'
        if val[2].endswith(pat):
            return key
    return None

# =============================================================================
# TASK 3.a: Mapping Benzene Ring Section
# =============================================================================
def map_benzene_ring_section(section: List[List[Any]],
                             final: List[str],
                             martini_dict: Dict[str, List[Any]],
                             full_mapping: List[List[List[Any]]]) -> List[str]:
    """
    Map a benzene ring section (section type == 2) using a tree‐like algorithm.
    Implements the following steps:
      1. Build array0: candidate atoms with type 'c/C' that have outer connections where the foreign atom's
         section type != 0, and that are unmapped.
      2. Pair atoms in array0 (using inner connections) and assign bead "TC5e" + random string.
      3. Build array1: atoms with no outer connections or outer connections of size 3+ or 2 that aren't CO
         that are unmapped;
     	if the amount of atoms in this array is odd: 
            first find any atom that is not connected to any other atoms in this array and remove it. If found then go straight to the pairing 
            if not found, find any atom that is only connecting to exactly 1 other atom in this array and remove it, if found then go to the pairing 
            if not found, throw an error saying benzene ring odd parity is bugged
            pair them (leftovers → array2) and assign beads ("TC5" if both C, "TN6a" if one N).
      4. Build array3: atoms that have outer connections where for at least one connection the foreign
         atom’s section is type 0 and of length 1.
         Then for each such atom:
            - If it has exactly one neighbor (via inner connections) found in array2, remove that neighbor
              from array2 and assign a bead based on the foreign atom type (using the first outer connection):
                  O: "SN6", N: "SN6d", S: "SC6", Cl: "SX3", I: "X1", C: "SC4"
            - If it has two neighbors, choose one (preferring a neighbor adjacent to an already mapped atom)
              and assign similarly.
            - If it has no neighbor in array2, use the bond order of the outer connection to decide:
                  if O with bond==2: "TN6a", if O with bond==1: "TN6", if C: "TC4" and more 
      5. Build array4: atoms with outer connections that are unmapped and for which the foreign section is
         of length 2; for each, if the foreign atom is O and its connected atom in that section is C,
         assign bead "SN2a" + random string.
      6. Build array5: every remaining unmapped atom; pair these using inner connections and assign beads
         ("TC5" if both C, "TN6a" if one N).
      7. Final check: ensure every atom is mapped.
    """
    # ---------------------- Step 1: Build array0 ----------------------
    array0 = []
    seen0 = set()
    for atom in section:
        idx = atom[0]
        # base criteria: carbon, unmapped, has any outer (foreign) connection
        if atom[1].lower() == 'c' and final[idx] == "" and atom[3]:
            valid_foreign = any(full_mapping[tup[0]][tup[1]][2] != 0 for tup in atom[3])
            if not valid_foreign:
                continue
    
            # look for an inner neighbor that also meets those criteria
            for inner_tup in atom[4]:
                nbr_global = section[inner_tup[0]][0]
                if nbr_global in seen0 or nbr_global == idx:
                    continue
    
                nbr_obj = next(a for a in section if a[0] == nbr_global)
                if (nbr_obj[1].lower() == 'c'
                    and final[nbr_global] == ""
                    and nbr_obj[3]
                    and any(full_mapping[t[0]][t[1]][2] != 0 for t in nbr_obj[3])):
                    # both carbons form a valid pair—add them once
                    array0.extend([idx, nbr_global])
                    seen0.update([idx, nbr_global])
                    break
    
    # array0 should now contain paired indices: length 2, 4, or 6
    
    # ---------------------- Step 2: Pair atoms in array0 ----------------------
    pairs0 = []
    used0 = set()
    for atom in section:
        if atom[0] in array0 and atom[0] not in used0:
            for tup in atom[4]:
                nbr_local = tup[0]
                if 0 <= nbr_local < len(section):
                    neighbor = section[nbr_local]
                    if neighbor[0] in array0 and neighbor[0] not in used0:
                        pairs0.append((atom[0], neighbor[0]))
                        used0.add(atom[0])
                        used0.add(neighbor[0])
                        break
    for p in pairs0:
        rstr = generate_random_string()
        # Select one of the atoms (here, the one with index p[0]) and get its first valid outer connection.
        atom_candidate = next(a for a in section if a[0] == p[0])
        bead_type = None
        for tup in atom_candidate[3]:
            foreign_atom = full_mapping[tup[0]][tup[1]]
            if foreign_atom[2] != 0:
                if foreign_atom[2] == 2:
                    bead_type = "TC5e"
                elif foreign_atom[2] == 1:
                    bead_type = "TC5"
                break
        # Default to "TC5e" if none found (should not happen if Step 1 is correct).
        if bead_type is None:
            bead_type = "TC5e"
        final[p[0]] = bead_type + rstr
        final[p[1]] = bead_type + rstr

    # ---------------------- Step 3: Build array1 ----------------------
    # Atoms that are unmapped and either have no outer connections OR have at least one outer connection
    # where the foreign section satisfies:
    #   - length >= 3, OR
    #   - length == 2 and the two atoms are not (one C and one O)
    def qualifies_foreign(tup):
        foreign_sec = full_mapping[tup[0]]
        if len(foreign_sec) >= 3:
            return True
        elif len(foreign_sec) == 2:
            # Get the uppercase types of both atoms in the foreign section.
            types = [foreign_sec[0][1].upper(), foreign_sec[1][1].upper()]
            if sorted(types) != ['C', 'O']:
                return True
        return False

    array1 = [atom[0] for atom in section 
              if final[atom[0]] == "" and 
              (len(atom[3]) == 0 or any(qualifies_foreign(tup) for tup in atom[3]))]
    # New parity correction: if array1 has an odd number of atoms.
    # Leftover atoms go to array2.
    array2 = []
    # try to remove any atom that is not connected to any other atom in array1.
    for idx in list(array1):  # iterate over a copy to avoid modification during iteration
        atom_obj = next(a for a in section if a[0] == idx)
        connection_count = sum(1 for tup in atom_obj[4] 
                               if 0 <= tup[0] < len(section) and section[tup[0]][0] in array1)
        if connection_count == 0:
            array2.append(idx)
            array1.remove(idx)
            removed = True
            break
    if len(array1) % 2 == 1:
        removed = False
        # First, try to remove any atom that is not connected to any other atom in array1.
        for idx in list(array1):  # iterate over a copy to avoid modification during iteration
            atom_obj = next(a for a in section if a[0] == idx)
            connection_count = sum(1 for tup in atom_obj[4] 
                                   if 0 <= tup[0] < len(section) and section[tup[0]][0] in array1)
            if connection_count == 0:
                array2.append(idx)
                array1.remove(idx)
                removed = True
                break
        # If none found, try to remove any atom that is connected to an inner atom 
        # that is unmapped and is not in array1.
        if not removed:
            for idx in list(array1):
                atom_obj = next(a for a in section if a[0] == idx)
                # Look for any inner connection to an unmapped atom not in array1.
                for tup in atom_obj[4]:
                    if 0 <= tup[0] < len(section):
                        neighbor = section[tup[0]]
                        # Check that the neighbor is unmapped and its global index is not in array1.
                        if final[neighbor[0]] == "" and neighbor[0] not in array1:
                            array2.append(idx)
                            array1.remove(idx)
                            removed = True
                            break
                if removed:
                    break
        if not removed:
            raise ValueError("Benzene ring odd parity is bugged!")
    # ---------------------- Step 4: Pair atoms from array1 ----------------------
    # Revised pairing for array1 according to rules described.
    pairs1 = []
    
    if len(array1) == 2:
        # If exactly two atoms remain, pair them.
        pairs1.append((array1[0], array1[1]))
    
    elif len(array1) == 4:
        # Try to find an atom that is connected to exactly 1 other atom in array1.
        candidate_pair = None
        for a_global in array1:
            a_atom = next(atom for atom in section if atom[0] == a_global)
            # Gather neighbors from inner connections, adjusting using the neighbor's local position.
            connected_neighbors = []
            for tup in a_atom[4]:
                nbr_local = tup[0]
                if 0 <= nbr_local < len(section):
                    neighbor = section[nbr_local]
                    if neighbor[0] in array1:
                        connected_neighbors.append(neighbor[0])
            if len(connected_neighbors) == 1:
                candidate_pair = (a_global, connected_neighbors[0])
                break
        if candidate_pair is not None:
            pairs1.append(candidate_pair)
            remaining = [x for x in array1 if x not in candidate_pair]
            # The remaining should be two atoms; pair them.
            if len(remaining) == 2:
                pairs1.append((remaining[0], remaining[1]))
            else:
                # If something unexpected happens, pair sequentially.
                for i in range(0, len(remaining), 2):
                    pairs1.append((remaining[i], remaining[i+1]))
        else:
            # Fallback: sort array1 by section order and pair sequentially.
            sorted_array1 = sorted(array1, key=lambda x: next(i for i, atom in enumerate(section) if atom[0] == x))
            for i in range(0, len(sorted_array1), 2):
                pairs1.append((sorted_array1[i], sorted_array1[i+1]))
    
    elif len(array1) == 6:
        # First, try to find any atom that connects to at least one other in array1.
        candidate_pair = None
        for a_global in array1:
            a_atom = next(atom for atom in section if atom[0] == a_global)
            connected_neighbors = []
            for tup in a_atom[4]:
                nbr_local = tup[0]
                if 0 <= nbr_local < len(section):
                    neighbor = section[nbr_local]
                    if neighbor[0] in array1:
                        connected_neighbors.append(neighbor[0])
            if connected_neighbors:
                candidate_pair = (a_global, connected_neighbors[0])
                break
        if candidate_pair is not None:
            pairs1.append(candidate_pair)
            remaining = [x for x in array1 if x not in candidate_pair]
            # Now, remaining should have 4 atoms. Process them as in the 4-case.
            candidate_pair_4 = None
            for a_global in remaining:
                a_atom = next(atom for atom in section if atom[0] == a_global)
                connected_neighbors = []
                for tup in a_atom[4]:
                    nbr_local = tup[0]
                    if 0 <= nbr_local < len(section):
                        neighbor = section[nbr_local]
                        if neighbor[0] in remaining:
                            connected_neighbors.append(neighbor[0])
                if len(connected_neighbors) == 1:
                    candidate_pair_4 = (a_global, connected_neighbors[0])
                    break
            if candidate_pair_4 is not None:
                pairs1.append(candidate_pair_4)
                remaining = [x for x in remaining if x not in candidate_pair_4]
                if len(remaining) == 2:
                    pairs1.append((remaining[0], remaining[1]))
                else:
                    for i in range(0, len(remaining), 2):
                        pairs1.append((remaining[i], remaining[i+1]))
            else:
                sorted_remaining = sorted(remaining, key=lambda x: next(i for i, atom in enumerate(section) if atom[0]==x))
                for i in range(0, len(sorted_remaining), 2):
                    pairs1.append((sorted_remaining[i], sorted_remaining[i+1]))
        else:
            # Fallback: sort array1 and pair sequentially.
            sorted_array1 = sorted(array1, key=lambda x: next(i for i, atom in enumerate(section) if atom[0] == x))
            for i in range(0, len(sorted_array1), 2):
                pairs1.append((sorted_array1[i], sorted_array1[i+1]))
    
    # pairs1 now contains the paired candidates following the specified rules.
    # Assign beads for each pair in pairs1.
    for pair in pairs1:
        rstr = generate_random_string()
        atom1 = next(a for a in section if a[0] == pair[0])
        atom2 = next(a for a in section if a[0] == pair[1])
        if atom1[1].lower() == 'c' and atom2[1].lower() == 'c':
            final[pair[0]] = "TC5" + rstr
            final[pair[1]] = "TC5" + rstr
        elif 'n' in (atom1[1].lower(), atom2[1].lower()):
            final[pair[0]] = "TN6a" + rstr
            final[pair[1]] = "TN6a" + rstr
        else:
            final[pair[0]] = "TC5" + rstr
            final[pair[1]] = "TC5" + rstr

    # ---------------------- Step 5: array2 remains from array1 pairing ----------------------
    # (array2 holds atoms that could not be paired from array1.)

    # ---------------------- Step 6: Build array3 ----------------------
    # For every atom with an outer connection (index3 non-empty) that is still unmapped,
    # and for at least one outer connection the foreign atom (in full_mapping) has section type 0 and the foreign section has length 1.
    array3 = []
    for atom in section:
        if len(atom[3]) != 0 and final[atom[0]] == "":
            for tup in atom[3]:
                foreign_sec = full_mapping[tup[0]]
                foreign_atom = foreign_sec[tup[1]]
                # require section type 0, length 1, and unmapped
                if foreign_atom[2] == 0 and len(foreign_sec) == 1 and final[foreign_atom[0]] == "":
                    array3.append(atom[0])
                    break

    # ---------------------- Step 7: Process array3 ----------------------
    # For every atom in array3, check inner connections for neighbors that satisfy:
    # either the neighbor is in array2 OR has an outer connection to a section that meets one of the following:
    # (a) foreign section length >= 3, or
    # (b) foreign section length == 2 and the two atoms in that section are not one C and one O.
    for a_idx in array3:
        atom = next(a for a in section if a[0] == a_idx)
        nbr_candidates = []
        for tup in atom[4]:
            nbr_local = tup[0]
            bond_order = tup[1]
            if 0 <= nbr_local < len(section):
                neighbor = section[nbr_local]
                qualifies = False
                # Condition 1: directly in array2.
                if neighbor[0] in array2:
                    qualifies = True
                else:
                    # Check neighbor's outer connections.
                    for outer in neighbor[3]:
                        foreign_sec = full_mapping[outer[0]]
                        if len(foreign_sec) >= 3:
                            qualifies = True
                            break
                        elif len(foreign_sec) == 2:
                            # Get the types of both atoms in the foreign section.
                            types = [a[1].upper() for a in foreign_sec]
                            # Only disqualify if they are one C and one O.
                            if sorted(types) != ['C', 'O']:
                                qualifies = True
                                break
                if qualifies and final[neighbor[0]] == "":
                    nbr_candidates.append((neighbor[0], bond_order))
        rstr = generate_random_string()
        foreign_tup = atom[3][0]  # Use the first outer connection.
        foreign_sec = full_mapping[foreign_tup[0]]
        foreign_atom = foreign_sec[foreign_tup[1]]
        if len(nbr_candidates) == 1:
            candidate = nbr_candidates[0][0]
            if candidate in array2:
                array2.remove(candidate)
            bond_order = atom[3][0][2]
            key = pick_bead_key(martini_dict,
                                foreign_atom[1].upper(),
                                bond_order,
                                kind='S')
            bead = (key + rstr) if key else ''
            final[atom[0]] = bead
            final[candidate] = bead
            final[foreign_atom[0]] = bead
            
        elif len(nbr_candidates) == 2:
            # If two candidates qualify, choose one—preferring a candidate adjacent to an already mapped atom.
            chosen = None
            for candidate, _ in nbr_candidates:
                neighbor_atom = next(a for a in section if a[0] == candidate)
                if any(final[section[tup[0]][0]] != "" for tup in neighbor_atom[4]):
                    chosen = candidate
                    break
            if chosen is None:
                chosen = nbr_candidates[0][0]
            if chosen in array2:
                array2.remove(chosen)
            bond_order = atom[3][0][2]
            key = pick_bead_key(martini_dict,
                                foreign_atom[1].upper(),
                                bond_order,
                                kind='S')
            bead = (key + rstr) if key else ''
            final[atom[0]] = bead
            final[chosen] = bead
            final[foreign_atom[0]] = bead
            
        else:
            # No qualifying neighbor: decide based on the bond order of the first outer connection.
            bond_order = atom[3][0][2]
            rstr = generate_random_string()
            key = pick_bead_key(martini_dict, foreign_atom[1].upper(), bond_order, kind='T')
            bead = (key + rstr) if key else ''
            final[atom[0]] = bead
            final[foreign_atom[0]] = bead
                
    # ---------------------- Step 8: Build array4 ----------------------
    # Atoms with outer connections still unmapped where for at least one outer connection
    # the foreign section has length 2.
    array4 = []
    for atom in section:
        if len(atom[3]) != 0 and final[atom[0]] == "":
            for tup in atom[3]:
                foreign_sec = full_mapping[tup[0]]
                foreign_atom = foreign_sec[tup[1]]
                # require section type 0, length 2, and unmapped foreign atom
                if foreign_atom[2] == 0 and len(foreign_sec) == 2 and final[foreign_atom[0]] == "":
                    array4.append(atom[0])
                    break
    
    # ---------------------- Step 9: Process array4 ----------------------
    # Only map SN2a if the number of remaining unmapped atoms in this section,
    # after subtracting the size of array4, is an even number.
    unmapped_count = len([atom for atom in section if final[atom[0]] == ""])
    if (unmapped_count - len(array4)) % 2 == 0:
        for a_idx in array4:
            atom = next(a for a in section if a[0] == a_idx)
            foreign_tup = atom[3][0]  # take the first outer connection
            foreign_sec = full_mapping[foreign_tup[0]]
            foreign_atom = foreign_sec[foreign_tup[1]]
            rstr = generate_random_string()
            if foreign_atom[1].upper() == 'O' and atom[1].upper() == 'C':
                final[atom[0]] = "SN2a" + rstr
                final[foreign_atom[0]] = "SN2a" + rstr
                if foreign_atom[4]:
                    nbr_local = foreign_atom[4][0][0]
                    nbr_in_foreign = foreign_sec[nbr_local]
                    final[nbr_in_foreign[0]] = "SN2a" + rstr
    else:
        # If the condition is not met, we skip mapping SN2a.
        pass

    # ---------------------- Step 10: Build array5 ----------------------
    # Every remaining unmapped atom.
    array5 = [atom[0] for atom in section if final[atom[0]] == ""]
    pairs5 = []
    # Case: six remaining
    if len(array5) == 6:
        used = set()
        for idx in array5:
            if idx in used:
                continue
            atom = next(a for a in section if a[0] == idx)
            for tup in atom[4]:  # inner connections
                nbr = section[tup[0]][0]
                if nbr in array5 and nbr not in used:
                    pairs5.append((idx, nbr))
                    used.add(idx)
                    used.add(nbr)
                    break
        # pair the last two
        rem = [i for i in array5 if i not in used]
        if len(rem) == 2:
            pairs5.append((rem[0], rem[1]))
    # Case: four remaining
    elif len(array5) == 4:
        used = set()
        # find a node with an already‑mapped neighbor
        for idx in array5:
            if idx in used:
                continue
            atom = next(a for a in section if a[0] == idx)
            for tup in atom[4]:
                nbr = section[tup[0]][0]
                if nbr in array5 and final[nbr] != "" and idx not in used and nbr not in used:
                    pairs5.append((idx, nbr))
                    used.add(idx)
                    used.add(nbr)
                    break
        # pair the remaining two
        rem = [i for i in array5 if i not in used]
        if len(rem) == 2:
            pairs5.append((rem[0], rem[1]))
    # Case: two remaining
    elif len(array5) == 2:
        i1, i2 = array5
        # if they’re directly connected, great:
        atom1 = next(a for a in section if a[0] == i1)
        if any(section[nbr_local][0] == i2 for nbr_local, _ in atom1[4]):
            pairs5 = [(i1, i2)]
        else:
            # fallback: reuse the first original pair from Step 2
            if not pairs0:
                # no original pair to borrow—just pair arbitrarily
                pairs5 = [(i1, i2)]
            else:
                orig_a, orig_b = pairs0[0]
    
                # 1) Unmap that bead entirely
                bead_tag = final[orig_a]
                for j, v in enumerate(final):
                    if v == bead_tag:
                        final[j] = ""
    
                # 2) Recompute who’s unmapped now (should be 4 atoms)
                unmapped = {a[0] for a in section if final[a[0]] == ""}
    
                # 3) For each of the two originals, pair it with one of its as-yet-unmapped neighbors
                for orig in (orig_a, orig_b):
                    atom = next(a for a in section if a[0] == orig)
                    for nbr_local, _ in atom[4]:
                        nbr_global = section[nbr_local][0]
                        if nbr_global in unmapped and nbr_global != orig:
                            pairs5.append((orig, nbr_global))
                            # mark both as used
                            unmapped.remove(orig)
                            unmapped.remove(nbr_global)
                            break
    # Map generated pairs
    for pair in pairs5:
        rstr = generate_random_string()
        atom1 = next(a for a in section if a[0] == pair[0])
        atom2 = next(a for a in section if a[0] == pair[1])
        if atom1[1].lower() == 'c' and atom2[1].lower() == 'c':
            final[pair[0]] = "TC5" + rstr
            final[pair[1]] = "TC5" + rstr
        elif 'n' in (atom1[1].lower(), atom2[1].lower()):
            final[pair[0]] = "TN6a" + rstr
            final[pair[1]] = "TN6a" + rstr

    # ---------------------- Final Check ----------------------
    for atom in section:
        if final[atom[0]] == "":
            raise ValueError("Benzene ring section not fully mappable!")
    return final


# =============================================================================
# TASK 3.b: Mapping Non-Benzene 6-Ring Section
# =============================================================================
def map_nonbenzene_6_ring_section(
    section: List[List[Any]],
    final: List[str],
    martini_dict: Dict[str, List[Any]],
    full_mapping: List[List[List[Any]]]
) -> List[str]:
    """
    Revised mapping function for non-benzene 6-ring sections incorporating outer
    (foreign) connections.

    The algorithm proceeds in four major steps:
    
    1. Process ring–ring border nodes: For each unmapped node having an outer connection
       to a group with section != 0, check its inner neighbors:
         - Case B1: If exactly one inner neighbor qualifies and one of the nodes
           has exactly 2 outer connections with one of the non-shared foreign connections
           (foreign section == 0 and length == 1), assign the candidate, its inner neighbor,
           and that foreign atom to a bead (using the provided rules snippet).
         - Case B2: If both nodes have 2 such foreign connections, raise an error.
         - Case B3: Otherwise, simply assign the two atoms a "TC5" bead.
    
    2. Process double–bond connections: For each unmapped pair connected by an inner bond of
       order 2, check their outer connections:
         - Case D1: If one atom has a qualifying foreign connection (section == 0, size == 1),
           assign all three atoms to a bead using the provided snippet.
         - Case D2: If both atoms have such foreign connections, map each with its foreign
           partner using alternate rules.
         - Case D3: Else assign the pair a "TC5" bead.
    
    3. Process remaining unmapped nodes that have an outer connection to a non–ring (section 0)
       group of size 1:
         - Build an auxiliary array from eligible inner neighbors (based on two steps ahead).
         - If no eligible neighbor is found, map the node and its chosen foreign neighbor as a pair.
         - If exactly one is found, map the node, that neighbor, and the foreign neighbor together.
         - If two are found, raise an error.
    
    4. Final fallback: Count remaining unmapped atoms.
         - If exactly 6 remain, use the original (oxygen-based) mapping.
         - Otherwise, for counts 5, 4, 3, 2, or 1, fixed
    """
    
    # Helper function: Get foreign atom from an outer connection tuple.
    def get_foreign_info(tup):
        foreign_sec = full_mapping[tup[0]]
        foreign_atom = foreign_sec[tup[1]]
        return foreign_atom
    
    # Helper function: Check if an outer connection qualifies as non-ring (section 0 with size 1)
    def qualifies_foreign(tup):
        foreign_sec = full_mapping[tup[0]]
        foreign_atom = foreign_sec[tup[1]]
        return (foreign_atom[2] == 0 and len(foreign_sec) == 1)
    
    # --- Step A: Process ring–ring border nodes
    for atom in section:
        if final[atom[0]] == "" and any(full_mapping[t[0]][t[1]][2] != 0 for t in atom[3]):
            candidate = atom
            matching_neighbors = []
            # Check candidate's inner neighbors for one that has an outer connection to the same ring
            for t in candidate[4]:
                nbr_local = t[0]
                if 0 <= nbr_local < len(section):
                    neighbor = section[nbr_local]
                    if final[neighbor[0]] == "":
                        # Does neighbor have an outer connection with the same section type as candidate?
                        for ft in neighbor[3]:
                            foreign_sec = full_mapping[ft[0]]
                            foreign_atom = foreign_sec[ft[1]]
                            if foreign_atom[2] == candidate[2]:
                                matching_neighbors.append(neighbor)
                                break
            if len(matching_neighbors) == 0:
                print("ring-ring neighbors connected through only 1 atom, is ignore")
            elif len(matching_neighbors) == 1:
                inner_neighbor = matching_neighbors[0]
                cand_qual = [t for t in candidate[3] if qualifies_foreign(t)]
                neigh_qual = [t for t in inner_neighbor[3] if qualifies_foreign(t)]
                if (len(cand_qual) == 1 and len(neigh_qual) == 1):
                    # Case B2: Too complex intermediate section
                    # copy from Case D2: Both nodes have a qualifying foreign connection.
                    # Use each node's first qualifying connection
                    fc1 = cand_qual[0]
                    fc2 = neigh_qual[0]
                    foreign_atom1 = get_foreign_info(fc1)
                    foreign_atom2 = get_foreign_info(fc2)
                    # Use the same bond order from one of the nodes
                    bond_order = atom[3][0][2]
                    key1 = pick_bead_key(martini_dict, foreign_atom1[1].upper(), bond_order, kind='T')
                    key2 = pick_bead_key(martini_dict, foreign_atom2[1].upper(), bond_order, kind='T')
                    rstr = generate_random_string()
                    bead1 = (key1 + rstr) if key1 else ''
                    rstr = generate_random_string()
                    bead2 = (key2 + rstr) if key2 else ''
                    final[candidate[0]] = bead1
                    final[foreign_atom1[0]] = bead1
                    final[inner_neighbor[0]] = bead2
                    final[foreign_atom2[0]] = bead2
                elif (len(cand_qual) == 1) or (len(neigh_qual) == 1):
                    # Case B1: Use the one (candidate or neighbor) that has 2 outer connections.
                    rstr = generate_random_string()
                    fc = cand_qual[0] if len(cand_qual) == 1 else neigh_qual[0]
                    foreign_atom = get_foreign_info(fc)
                    bond_order = candidate[3][0][2]  # using candidate's first outer connection's bond order
                    key = pick_bead_key(martini_dict,
                                        foreign_atom[1].upper(),
                                        bond_order,
                                        kind='S')
                    bead = (key + rstr) if key else ''
                    final[candidate[0]] = bead
                    final[inner_neighbor[0]] = bead
                    final[foreign_atom[0]] = bead
                else:
                    # Case B3: Assign the two atoms as TC5.
                    rstr = generate_random_string()
                    bead = "TC5" + rstr
                    final[candidate[0]] = bead
                    final[inner_neighbor[0]] = bead
            else:
                # Ambiguous: more than one matching inner neighbor → skip this atom
                print(f"Ambiguous ring–ring border for atom {atom[0]}; skipping")
                continue
    # --- Step B: Process double–bond connections (inner bonds with order 2) ---
    processed_double = set()
    for atom in section:
        if final[atom[0]] != "":
            continue
    
        for t in atom[4]:
            nbr_local, bo = t
            if bo != 2 and bo != 1.5:
                continue
    
            if not (0 <= nbr_local < len(section)):
                continue
    
            neighbor = section[nbr_local]
            a_idx = atom[0]
            b_idx = neighbor[0]
    
            # skip if already mapped or already processed
            if final[b_idx] != "" or (a_idx, b_idx) in processed_double or (b_idx, a_idx) in processed_double:
                continue
    
            # pull the bond order directly from the inner‐bond tuple
            bond_order = bo
    
            # find qualifying foreign connections on each
            cand_qual  = [x for x in atom[3]     if qualifies_foreign(x)]
            neigh_qual = [x for x in neighbor[3] if qualifies_foreign(x)]
    
            # Case D1: exactly one side qualifies
            if len(cand_qual) == 1 or len(neigh_qual) == 1:
                fc = cand_qual[0] if len(cand_qual) == 1 else neigh_qual[0]
                foreign = get_foreign_info(fc)
                key = pick_bead_key(martini_dict,
                                    foreign[1].upper(),
                                    bond_order,
                                    kind='S')
                bead = (key + generate_random_string()) if key else ''
                final[a_idx] = bead
                final[b_idx] = bead
                final[foreign[0]] = bead
    
            # Case D2: both sides qualify
            elif len(cand_qual) == 1 and len(neigh_qual) == 1:
                fa = get_foreign_info(cand_qual[0])
                fb = get_foreign_info(neigh_qual[0])
                key1 = pick_bead_key(martini_dict, fa[1].upper(), bond_order, kind='T')
                key2 = pick_bead_key(martini_dict, fb[1].upper(), bond_order, kind='T')
                bead1 = (key1 + generate_random_string()) if key1 else ''
                bead2 = (key2 + generate_random_string()) if key2 else ''
                final[a_idx] = bead1
                final[fa[0]] = bead1
                final[b_idx] = bead2
                final[fb[0]] = bead2
    
            # Case D3: neither side qualifies – choose bead type based on atom types
            else:
                # get uppercase element symbols for each of the two atoms
                elem_a = atom[1].upper()
                elem_b = neighbor[1].upper()
    
                # if both are C → TC5
                if elem_a == 'C' and elem_b == 'C':
                    bead = "TC5" + generate_random_string()
    
                # if one is C and the other is N → TN6a
                elif (elem_a == 'C' and elem_b == 'N') or (elem_a == 'N' and elem_b == 'C'):
                    bead = "TN6a" + generate_random_string()
    
                # otherwise, fall back to TC5
                else:
                    bead = "TC5" + generate_random_string()
    
                final[a_idx] = bead
                final[b_idx] = bead
    
            processed_double.add((a_idx, b_idx))
            # once we’ve assigned for this atom, break out to next atom
            break
    # --- Step C: Process remaining unmapped nodes with 1 or 2 non-ring size-1 foreign connections ---
    for i in range(2):
        for atom in section:
            if final[atom[0]] == "":
                # collect all size-1 non-ring outer tuples whose foreign atom is still unmapped
                f_tups = [t for t in atom[3]
                          if qualifies_foreign(t) and final[full_mapping[t[0]][t[1]][0]] == ""]
                # case: two size-1 foreigns → check all-C and assign SC3
                if len(f_tups) == 2:
                    # unpack the two foreign‐tuples (sec_idx, loc_idx, bond_order)
                    ft1, ft2 = f_tups
                    fa1 = get_foreign_info(ft1)  # [global_idx, element, ...]
                    fa2 = get_foreign_info(ft2)
                    elem = atom[1].upper()
                    e1, e2 = fa1[1].upper(), fa2[1].upper()
                
                    # helper to pick and tag a bead
                    def assign(bead_key: str):
                        tag = bead_key + generate_random_string()
                        final[atom[0]] = tag
                        final[fa1[0]]  = tag
                        final[fa2[0]]  = tag
                
                    if elem == 'C':
                        # case: one foreign is C, the other is something else
                        if e1 == 'C' or e2 == 'C':
                            # pick the non‑C partner
                            chosen_ft = ft1 if e1 != 'C' else ft2
                            fa        = get_foreign_info(chosen_ft)
                            bond      = chosen_ft[2]  # use bond_order=1 here per your spec
                            key       = pick_bead_key(martini_dict, fa[1].upper(), bond, kind='S')
                            if not key:
                                raise ValueError(f"No S‑key for {fa[1]}")
                            assign(key)
                
                        # both O → SN5a
                        elif e1 == 'O' and e2 == 'O':
                            assign('SN5a')
                
                        # both N → SN1
                        elif e1 == 'N' and e2 == 'N':
                            assign('SN1')
                        # both F → SX4e
                        elif e1 == 'F' and e2 == 'F':
                            assign('SX4e')
                
                        else:
                            raise ValueError("2 size‑1 neighbors incompatible for central C")
                
                    elif elem == 'O':
                        # only C–C allowed → SN3a
                        if e1 == 'C' and e2 == 'C':
                            assign('SN3a')
                        else:
                            raise ValueError("2 size‑1 neighbors incompatible for central O")
                
                    else:
                        raise ValueError("2 size‑1 neighbors but central atom not C/O")
                
                    continue  # move on after handling this case
                # case: single size-1 foreign → proceed with original single-tuple logic
                if len(f_tups) == 1:
                    primary = f_tups[0]
                    foreign_atom = get_foreign_info(primary)
                    bond_order = primary[2]
                    # build array1 of candidate inner neighbors
                    array1 = []
                    for t in atom[4]:
                        nbr = section[t[0]]
                        if final[nbr[0]] == "":
                            # check if all of neighbor's inner connections remain unmapped
                            inner_unmapped = all(
                                final[section[t2[0]][0]] == "" for t2 in nbr[4]
                            )
                            # gather neighbor's size-1 non-ring foreign tuples unmapped
                            nf_tups = [x for x in nbr[3]
                                       if qualifies_foreign(x) and final[full_mapping[x[0]][x[1]][0]] == ""]
                            # include if neighbor is C and not inner_unmapped, with tuple rules
                            if nbr[1].upper() == 'C' and not inner_unmapped:
                                if len(nf_tups) == 1 or (len(nbr[3]) == 2 and len(nf_tups) == 0) or len(nf_tups) == 0:
                                    array1.append(nbr)
                    # map based on how many neighbors found
                    if len(array1) == 0:
                        rstr = generate_random_string()
                        key = pick_bead_key(martini_dict, foreign_atom[1].upper(), bond_order, kind='T')
                        bead = (key + rstr) if key else ''
                        final[atom[0]] = bead
                        final[foreign_atom[0]] = bead
                    elif len(array1) == 1:
                        rstr = generate_random_string()
                        key = pick_bead_key(martini_dict,
                                            foreign_atom[1].upper(),
                                            bond_order,
                                            kind='S')
                        bead = (key + rstr) if key else ''
                        final[atom[0]] = bead
                        final[array1[0][0]] = bead
                        final[foreign_atom[0]] = bead
                    elif len(array1) == 2:
                        if i == 1:
                            raise ValueError("Non–ring section mapping too complex (foreign/inner connection).")
    # --- Step C.5: Special catch‐all for any unmapped inner atom that is not C or O,
    #                 but has a single‐sized non‑ring (section 0) foreign neighbor which is C.
    #                 Map both atoms to a T‐type bead (bond order = 1).
    for atom in section:
        if final[atom[0]] != "":
            continue
        # If atom is neither C nor O, attempt to see if it has a qualifying outer neighbor that is C
        elem = atom[1].upper()
        if elem in ('C', 'O'):
            continue
    
        # Look for a non‑ring foreign connection of size 1 (qualifies_foreign) where the foreign atom is C
        for ft in atom[3]:
            # ft = (foreign_section_index, foreign_atom_local_index, bond_order)
            foreign_atom = get_foreign_info(ft)
            if qualifies_foreign(ft) and foreign_atom[1].upper() == 'C':
                # We have an inner atom (elem ≠ C/O) and a foreign carbon
                # Use bond_order = 1 and kind='T'
                key = pick_bead_key(martini_dict, elem, 1, kind='T')
                bead = (key + generate_random_string()) if key else ''
                # Assign both inner atom and foreign carbon to the same bead
                final[atom[0]] = bead
                final[foreign_atom[0]] = bead
                break   # move on to next atom
                
    # --- Step D: Final fallback mapping for remaining unmapped atoms
    remaining = [atom for atom in section if final[atom[0]] == ""]
    count = len(remaining)
    # ensure only C or O remain
    if count in (5,4,3,2,1):
        if any(atom[1].upper() not in ('C','O') for atom in remaining):
            raise ValueError("non-benzene 6-ring section has non C/O atom")
    if count == 6:
        # Fall back to the original oxygen-based mapping.
        array_O = [atom[0] for atom in section if atom[1].upper() == "O" and final[atom[0]] == ""]
        if len(array_O) == 1:
            target_idx = array_O[0]
            target = next(a for a in section if a[0] == target_idx)
            a_str = generate_random_string()
            bead = "SN4a" + a_str
            final[target[0]] = bead
            for tup in target[4]:
                neighbor_idx = section[tup[0]][0]
                final[neighbor_idx] = bead
        elif len(array_O) == 2:
            for idx in array_O:
                a_str = generate_random_string()
                bead = "SN3a" + a_str
                final[idx] = bead
                atom_obj = next(a for a in section if a[0] == idx)
                for tup in atom_obj[4]:
                    neighbor_idx = section[tup[0]][0]
                    final[neighbor_idx] = bead
        array1 = [atom[0] for atom in section if final[atom[0]] == ""]
        if len(array1) == 3:
            a_str = generate_random_string()
            bead = "SC3" + a_str
            for idx in array1:
                final[idx] = bead
        elif len(array1) == 6:
            group1 = array1[:3]
            group2 = array1[3:]
            a_str1 = generate_random_string()
            bead1 = "SC3" + a_str1
            for idx in group1:
                final[idx] = bead1
            a_str2 = generate_random_string()
            bead2 = "SC3" + a_str2
            for idx in group2:
                final[idx] = bead2
        else:
            raise ValueError("Fallback mapping for non-benzene 6-ring with count {} not implemented (TODO).".format(count))
    elif count == 5:
        # find the C-atom with exactly one unmapped C neighbor
        cand = None
        for atom in remaining:
            if final[atom[0]] != "":
                continue
            nbrs = [section[t[0]] for t in atom[4]
                    if final[ section[t[0]][0] ] == ""]
            if len(nbrs) == 1 and nbrs[0][1].upper() == 'C':
                cand = (atom, nbrs[0]); break
        if not cand:
            raise ValueError("Non-benzene 6-ring: no singleton C neighbor for count=5")
        atom, nbr = cand
        old_bead = final[nbr[0]]
        # find all indices in final with old_bead
        targets = [i for i,v in enumerate(final) if v == old_bead]
        new_bead = old_bead
        if new_bead.startswith('T'):
            new_bead = 'S' + new_bead[1:]
        for i in targets:
            final[i] = new_bead
        final[atom[0]] = new_bead
        # now fall through to count=4 logic
        count = 4
        remaining = [atom for atom in section if final[atom[0]] == ""]
    if count == 4:
        # for each atom with only one unmapped neighbor, build type string
        for atom in remaining:
            idx = atom[0]
            if final[idx] != "":
                continue
            nbrs = [section[t[0]] for t in atom[4] if final[ section[t[0]][0] ] == ""]
            if len(nbrs) == 1:
                types = atom[1].upper() + nbrs[0][1].upper()
                # find T-keys with val[0]==3 or 5 matching types or reversed
                for key,val in martini_dict.items():
                    if key.startswith('T') and val[0] in (3,5) and (val[2] == types or val[2] == types[::-1]):
                        bead = key + generate_random_string()
                        final[atom[0]] = bead
                        final[nbrs[0][0]] = bead
                        break
    elif count == 3:
        # center atom has two unmapped neighbors
        center = next(atom for atom in remaining
                      if len([t for t in atom[4] if final[ section[t[0]][0] ] == ""]) == 2)
        nbrs = [section[t[0]] for t in center[4] if final[ section[t[0]][0] ] == ""]
        s = nbrs[0][1].upper() + center[1].upper() + nbrs[1][1].upper()
        for key,val in martini_dict.items():
            if key.startswith('S') and val[0] in (3,5) and (val[2] == s or val[2] == s[::-1]):
                bead = key + generate_random_string()
                final[center[0]] = bead
                for n in nbrs: final[n[0]] = bead
                break
        # if still unmapped, do the “merge” logic
        if final[center[0]] == "":
            merged = False

            # 1) look for a C-neighbor with a T-bead neighbor
            for nbr in nbrs:
                if nbr[1].upper() == 'C':
                    # gather all neighbor indices (inner + outer)
                    inner_idxs = [section[t[0]][0] for t in nbr[4]]
                    outer_idxs = [full_mapping[t[0]][t[1]][0] for t in nbr[3]]
                    all_idxs = inner_idxs + outer_idxs

                    t_beads = [final[i] for i in all_idxs if final[i].startswith('T')]
                    if t_beads:
                        old = t_beads[0]
                        new = 'S' + old[1:]
                        # reassign everywhere
                        for j, v in enumerate(final):
                            if v == old:
                                final[j] = new
                        final[nbr[0]] = new
                        merged = True
                        break

            # 2) if still not merged, fall back to any S-bead neighbor
            if not merged:
                for nbr in nbrs:
                    inner_idxs = [section[t[0]][0] for t in nbr[4]]
                    outer_idxs = [full_mapping[t[0]][t[1]][0] for t in nbr[3]]
                    all_idxs = inner_idxs + outer_idxs

                    s_beads = [final[i] for i in all_idxs if final[i].startswith('S')]
                    if s_beads:
                        bead = s_beads[0]
                        final[nbr[0]] = bead
                        merged = True
                        break

            # 3) if we still failed, error out
            if not merged:
                raise ValueError(
                    "Merging is not possible; mapping failed for count=3"
                )
            count = 2
            remaining = [atom for atom in section if final[atom[0]] == ""]
    if count == 2:
        a1,a2 = remaining
        types = a1[1].upper() + a2[1].upper()
        for key,val in martini_dict.items():
            if key.startswith('T') and val[0] in (3,5) and (val[2] == types or val[2] == types[::-1]):
                bead = key + generate_random_string()
                final[a1[0]] = bead
                final[a2[0]] = bead
                break
    elif count == 1:
        lone = remaining[0]
        if lone[1].upper() == 'O':
            raise ValueError("non benzene 6-ring has a lone Oxygen left")
        # else C
        # try find T-bead neighbor
        # check its mapped neighbors
        nbrs = [section[t[0]] for t in lone[4]]
        t_beads = [final[n[0]] for n in nbrs if final[n[0]].startswith('T')]
        if t_beads:
            old = t_beads[0]
            new = old
            if not new.startswith('S'): new = 'S' + new[1:]
            for i,v in enumerate(final):
                if v == old: final[i] = new
            final[lone[0]] = new
        else:
            s_beads = [final[n[0]] for n in nbrs if final[n[0]].startswith('S')]
            if s_beads:
                final[lone[0]] = s_beads[0]
            else:
                raise ValueError("non benzene 6-ring: cannot assign lone C bead")
    return final

# =============================================================================
# TASK 3.c: Mapping Non-Benzene 5-Ring Section
# =============================================================================
def map_nonbenzene_5_ring_section(section: List[List[Any]],
                                  final: List[str],
                                  martini_dict: Dict[str, List[Any]],
                                  full_mapping: List[List[List[Any]]]) -> List[str]:
    """
    Map a non-benzene 5-membered ring section.
    
    Pseudocode:
      (A) Build array0:
          For every atom that has:
            - element (index 1) is 'C' (or 'c')
            - final[atom[0]] is ""
            - AND it has at least one inner neighbor (from its inner connections at index 4)
              such that the neighbor's element is 'C' or 'N' AND the bond equals 2,
          add the pair (atom[0], neighbor[0]) as a tuple to array0.
      For each tuple in array0:
          - If both atoms are C: 
              if both atoms each have a foreign neighbor with the foreign sections of size 1: 
                  return an error saying 5-ring section is too difficult to be mapped and add a todo comment for later 
              if one atom has a foreign neighbor with the foreign section of size 1: 
                  assign a bead based on the foreign atom type (using the first outer connection): 
                      O: "SN6", N: "SN6d", S: "SC6", Cl: "SX3", I: "X1", C: "SC4", Br: SX2 
                      assign it + random string to all 3 atoms  
              else:
                assign bead "TC5" + random string. assign bead "TC5" + random string.
          - If one atom is N: assign bead "TN6a" + random string.
      
      (B) Build array1:
          For every atom that is unmapped (final[atom[0]] == "") and its element is not 'C', add its global index.
          Then:
           - If length of array1 is 1:
                * For the single atom:
                    - If element is S:
                        * If both of its first two inner neighbors are unmapped, generate a random string and assign bead "SC6" + that string to the atom and both neighbors.
                        * Else, assign "TC6" + random string only to the atom.
                    - If element is NH: 
                        * If both of its first two inner neighbors are unmapped, generate a random string and assign bead "TN6d" + that string to the atom and only 1 neighbor.
                        * Else, assign "TN6d" + random string to the only atom.
                    - If element is O: assign "TN4a" + random string to the atom and (only) to its first inner neighbor.
                    - If element is N and it has a foreign connection (len(atom[3]) > 0): 
                        Use the first outer connection to look up the foreign atom in full_mapping and assign "TN1" + random string to both.
           - If length of array1 is 2:
                * If both atoms have element O:
                    - Find a common neighbor (via inner connections) and assign "SN5a" + random string to both atoms and that common neighbor.
      
      (C) Build array2:
          For every atom that is unmapped and has element 'C', add its global index.
          We assume array2 will have either 5, 3, or 2 elements.
           - If length is 5:
                * Loop through array2, and find an atom that has:
                      - element 'C' (or 'c')
                      - at least one outer connection (len(atom[3]) != 0)
                  If found:
                      - If the foreign neighbor (from the first outer connection) has element 'C' and the size of its section is 1,
                        assign "SC3" + random string to that atom, the foreign neighbor, and its first inner neighbor.
                      - Else if the foreign neighbor has element 'O' and its section size is 1,
                        assign "SN6" + random string similarly.
                      - Else if the foreign neighbor’s section size is > 1, assign "TC3" + random string to that atom and its first inner neighbor.
                  If not found:
                      - take any connecting 2 atoms and assign TC3 + random string to those two atoms
                * After this, loop through array2 and remove atoms that are now mapped.
           - If the remaining size is 2 and all atoms have element 'C', assign "TC3" + random string.
           - If the remaining size is 3 and all atoms have element 'C', assign "SC3" + random string.
      
      (D) Final check: if any atom in the section remains unmapped, raise an error.
      
    Returns final.
    """
    # -----------------------
    # (A) Process double bonds among C and N atoms.
    candidate_indices: List[Tuple[int, int, str, Optional[bool]]] = []
    for pos, atom in enumerate(section):
        idx = atom[0]
        if final[idx] != "":
            continue
        elem = atom[1].upper()
        if elem in ('C', 'N'):
            for tup in atom[4]:  # inner connections
                bond = tup[1]
                if bond in (2, 1.5):
                    if elem == 'C':
                        candidate_indices.append((pos, idx, 'c', None))
                    else:
                        has_foreign = len(atom[3]) > 0
                        candidate_indices.append((pos, idx, 'n', has_foreign))
                    break
    # Post-process nitrogen candidates: if exactly 2 candidates of type 'n' exist,
    # remove the one that has a foreign neighbor if the other does not.
    candidate_n = [ci for ci in candidate_indices if ci[2] == 'n']
    if len(candidate_n) == 2:
        if candidate_n[0][3] and not candidate_n[1][3]:
            candidate_indices = [ci for ci in candidate_indices if not (ci[2] == 'n' and ci[1] == candidate_n[0][1])]
        elif candidate_n[1][3] and not candidate_n[0][3]:
            candidate_indices = [ci for ci in candidate_indices if not (ci[2] == 'n' and ci[1] == candidate_n[1][1])]
        else:
            # Define a helper to count the number of inner connections
            # from the candidate to any other candidate in candidate_n.
            def count_connections(candidate, candidate_list, section):
                count = 0
                global_idx = candidate[1]
                atom = next(a for a in section if a[0] == global_idx)
                candidate_globals = [c[1] for c in candidate_list]
                for tup in atom[4]:
                    nbr_local = tup[0]
                    if 0 <= nbr_local < len(section):
                        neighbor = section[nbr_local]
                        if neighbor[0] in candidate_globals:
                            count += 1
                return count
            count0 = count_connections(candidate_n[0], candidate_n, section)
            count1 = count_connections(candidate_n[1], candidate_n, section)
            if count0 == 1 and count1 != 1:
                candidate_indices = [ci for ci in candidate_indices if not (ci[2]=='n' and ci[1] == candidate_n[0][1])]
            else:
                candidate_indices = [ci for ci in candidate_indices if not (ci[2]=='n' and ci[1] == candidate_n[1][1])]
    elif len(candidate_n) == 1 and len(candidate_indices) % 2 == 1:
        # If there is only 1 'n' candidate and 0 2 or 4 'c' candidates (odd), remove the single 'n'.
        candidate_indices = [ci for ci in candidate_indices if ci[2] != 'n']
    # Sort candidates by their position in the section.
    candidate_indices.sort(key=lambda x: x[0])
    # Extract the ordered global indices.
    ordered_candidates = [ci[1] for ci in candidate_indices]
    
    # Sort candidates by their position in the section.
    candidate_indices.sort(key=lambda x: x[0])
    # Extract the ordered global indices.
    ordered_candidates = [ci[1] for ci in candidate_indices]
    def generate_perfect_matchings(candidates: List[int]) -> List[List[Tuple[int, int]]]:
        """Recursively generate all perfect matchings for the candidate list."""
        if not candidates:
            return [[]]
        matchings = []
        first = candidates[0]
        for i in range(1, len(candidates)):
            pair = (first, candidates[i])
            rest = candidates[1:i] + candidates[i+1:]
            for submatching in generate_perfect_matchings(rest):
                matchings.append([pair] + submatching)
        return matchings
    
    all_matchings = generate_perfect_matchings(ordered_candidates)
    
    # Filter valid matchings: each pair must be connected.
    valid_matchings = []
    for matching in all_matchings:
        is_valid = True
        for idx1, idx2 in matching:
            if final[idx1] != "" or final[idx2] != "":
                is_valid = False; 
                break
            # Look up atom1 using its global index.
            atom1 = next(a for a in section if a[0] == idx1)
            # Check if any inner connection from atom1 points to an atom whose global index equals idx2.
            if not any(section[tup[0]][0] == idx2 for tup in atom1[4]):
                is_valid = False
                break
        if is_valid:
            valid_matchings.append(matching)
    
    if valid_matchings:
        pairs = valid_matchings[0]
    else:
        pairs = []
    
    # Now, for each valid pair, assign the bead.
    for pair in pairs:
        rstr = generate_random_string()
        atom1 = next(a for a in section if a[0] == pair[0])
        atom2 = next(a for a in section if a[0] == pair[1])
        if atom1[1].upper() == "C" and atom2[1].upper() == "C":
            # Check for foreign neighbors with a foreign section of size 1.
            has_foreign1 = any(len(full_mapping[tup[0]]) == 1 for tup in atom1[3])
            has_foreign2 = any(len(full_mapping[tup[0]]) == 1 for tup in atom2[3])
            if has_foreign1 and has_foreign2:
                # TODO: Handle the case where both atoms have a qualifying foreign neighbor.
                raise ValueError("5-ring section is too difficult to be mapped! TODO: add proper handling for both foreign neighbors case.")
            elif has_foreign1 or has_foreign2:
                # One of the atoms has a qualifying foreign neighbor.
                if has_foreign1:
                    foreign_tup = next(tup for tup in atom1[3] if len(full_mapping[tup[0]]) == 1)
                else:
                    foreign_tup = next(tup for tup in atom2[3] if len(full_mapping[tup[0]]) == 1)
                foreign_sec = full_mapping[foreign_tup[0]]
                foreign_atom = foreign_sec[foreign_tup[1]]
                foreign_element = foreign_atom[1].upper()
                # Search martini_dict for keys with section value 2 and a representation that starts with "CC(".
                # Then extract the substring in the parentheses and compare it to foreign_element.
                chosen_key = None
                for key, value in martini_dict.items():
                    if value[0] == 2 and value[2].startswith("CC("):
                        open_index = value[2].find("(")
                        close_index = value[2].find(")")
                        if open_index != -1 and close_index != -1:
                            element_in_parenthesis = value[2][open_index+1:close_index].upper()
                            if element_in_parenthesis == foreign_element:
                                chosen_key = key
                                break
                if chosen_key is None:
                    chosen_key = "TC5"  # fallback default
                bead = chosen_key + rstr
                # Assign the chosen bead to both atoms in the pair and to the foreign atom.
                final[pair[0]] = bead
                final[pair[1]] = bead
                final[foreign_atom[0]] = bead
            else:
                # Default case: assign "TC5" bead with a random string.
                final[pair[0]] = "TC5" + rstr
                final[pair[1]] = "TC5" + rstr
        elif "N" in (atom1[1].upper(), atom2[1].upper()):
            # If one atom is nitrogen, assign "TN6a" bead with a random string.
            final[pair[0]] = "TN6a" + rstr
            final[pair[1]] = "TN6a" + rstr
    # -----------------------
    # (B) Build array1 for atoms not of type C.
    array1 = [atom[0] for atom in section if final[atom[0]] == "" and atom[1].upper() != "C"]
    if len(array1) == 1:
        rstr = generate_random_string()
        idx = array1[0]
        atom = next(a for a in section if a[0] == idx)
        if atom[1].upper() == "S":
            # Check inner neighbors: if both first two are unmapped, assign SC6 to all; else, assign TC6 only to atom.
            if len(atom[4]) >= 2:
                inn1 = section[atom[4][0][0]]
                inn2 = section[atom[4][1][0]]
                if final[inn1[0]] == "" and final[inn2[0]] == "":
                    bead = "SC6" + rstr
                    final[atom[0]] = bead
                    final[inn1[0]] = bead
                    final[inn2[0]] = bead
                else:
                    final[atom[0]] = "TC6" + rstr
            else:
                final[atom[0]] = "TC6" + rstr
        elif atom[1].upper() == "N" and len(atom[3]) > 0:
            foreign_tup = atom[3][0]
            foreign_atom = full_mapping[foreign_tup[0]][foreign_tup[1]]
            if final[foreign_atom[0]] == "":
                # Use the global mapping for foreign connection.
                rstr = generate_random_string()
                bead = "TN1" + rstr
                final[atom[0]] = bead
                # Lookup in full_mapping.
                final[foreign_atom[0]] = bead
        elif atom[1].upper() == "NH" or atom[1].upper() == "N":
            # Extra rule for NH:
            # If there are at least two inner connections and both of the first two inner neighbors are unmapped,
            # assign the bead "TN6d" + random string to the NH atom and (arbitrarily) to the first neighbor.
            if len(atom[4]) >= 2:
                inn1 = section[atom[4][0][0]]
                inn2 = section[atom[4][1][0]]
                if final[inn1[0]] == "" and final[inn2[0]] == "":
                    bead = "TN6d" + rstr
                    final[atom[0]] = bead
                    final[inn1[0]] = bead   # assign to only one neighbor (choose the first)
                else:
                    final[atom[0]] = "TN6d" + rstr
            else:
                final[atom[0]] = "TN6d" + rstr
        elif atom[1].upper() == "O":
            # Check if there are at least two inner neighbors.
            if len(atom[4]) >= 2:
                # Get global indices for the first two inner neighbors.
                left_idx = section[atom[4][0][0]][0]
                right_idx = section[atom[4][1][0]][0]
                # If both inner neighbors are unmapped, assign TN4a.
                if final[left_idx] == "" and final[right_idx] == "":
                    a_str = generate_random_string()
                    bead = "TN4a" + a_str
                    final[atom[0]] = bead
                    # According to your pseudocode, assign only the left neighbor.
                    final[left_idx] = bead
                else:
                    # Otherwise, assign TN2a to the atom only.
                    a_str = generate_random_string()
                    bead = "TN2a" + a_str
                    final[atom[0]] = bead
            else:
                # If there are fewer than 2 inner neighbors, default to TN2a.
                a_str = generate_random_string()
                bead = "TN2a" + a_str
                final[atom[0]] = bead
    elif len(array1) == 2:
        # If both atoms are O.
        atom1 = next(a for a in section if a[0] == array1[0])
        atom2 = next(a for a in section if a[0] == array1[1])
        if atom1[1].upper() == "O" and atom2[1].upper() == "O":
            common = None
            for tup in atom1[4]:
                n1 = section[tup[0]]
                for tup2 in atom2[4]:
                    n2 = section[tup2[0]]
                    if n1[0] == n2[0]:
                        common = n1[0]
                        break
                if common is not None:
                    break
            if common is not None:
                rstr = generate_random_string()
                final[atom1[0]] = "SN5a" + rstr
                final[atom2[0]] = "SN5a" + rstr
                final[common] = "SN5a" + rstr

    # -----------------------
    # (C) Build array2: for unmapped Cs.
    array2 = [atom[0] for atom in section if final[atom[0]] == "" and atom[1].upper() == "C"]
    
    if len(array2) == 5:
        # Try to handle the “one C with an outer” case exactly once.
        assigned = False
        for idx in array2:
            atom = next(a for a in section if a[0] == idx)
            if atom[3]:  # has at least one outer connection
                foreign_tup  = atom[3][0]
                foreign_sec  = full_mapping[foreign_tup[0]]
                foreign_atom = foreign_sec[foreign_tup[1]]
                rstr = generate_random_string()
    
                if foreign_atom[1].upper() == "C" and len(foreign_sec) == 1:
                    bead = "SC3" + rstr
                elif foreign_atom[1].upper() == "O" and len(foreign_sec) == 1:
                    bead = "SN6" + rstr
                else:
                    bead = "TC3" + rstr
    
                # assign to the trio: this C, its foreign neighbour, and its first inner neighbour
                final[atom[0]] = bead
                final[foreign_atom[0]] = bead
                if atom[4]:
                    inn = section[atom[4][0][0]]
                    final[inn[0]] = bead
    
                assigned = True
                break  # << don’t assign again!
    
        if not assigned:
            # fallback: pair any two connected Cs once, with TC3
            for i, idx1 in enumerate(array2):
                atom1 = next(a for a in section if a[0] == idx1)
                for idx2 in array2[i+1:]:
                    if any(section[t][0] == idx2 for t,_ in atom1[4]):
                        bead = "TC3" + generate_random_string()
                        final[idx1] = bead
                        final[idx2] = bead
                        assigned = True
                        break
                if assigned:
                    break
            if not assigned:
                raise ValueError(
                   "Could not find any two connecting atoms in 5‑ring array2 for TC3 assignment."
                )
    
        # now remove those atoms from array2
        array2 = [idx for idx in array2 if final[idx] == ""]
    
    # At this point, array2 contains only the *unmapped* Cs after above.
    # If exactly two remain, pair them.  If three remain, do SC3 on all three.
    if len(array2) == 2:
        bead = "TC3" + generate_random_string()
        for idx in array2:
            final[idx] = bead
        array2.clear()
    elif len(array2) == 3:
        bead = "SC3" + generate_random_string()
        for idx in array2:
            final[idx] = bead
        array2.clear()

    # -----------------------
    # (D) Final check with merging fallback
    for atom in section:
        idx = atom[0]
        if final[idx]=='' and atom[1].upper() == "C":
            # gather neighbors: inner + outer
            neigh = [section[t[0]][0] for t in atom[4]] + [full_mapping[t[0]][t[1]][0] for t in atom[3]]
            # try TC5 merge
            merged=False
            for prefix in ('TC5',):
                for n in neigh:
                    bead_str = final[n]
                    if bead_str.startswith(prefix):
                        # clear old
                        old = bead_str
                        ids = [i for i,v in enumerate(final) if v==old]
                        for i in ids: final[i] = ''
                        # pick new S-key
                        bond = atom[4][0][1] if atom[4] else None
                        key = pick_bead_key(martini_dict, atom[1].upper(), bond, kind='S')
                        if not key:
                            raise ValueError("Fallback mapping failed: no S-key for element")
                        bead = key + generate_random_string()
                        for i in ids: final[i]=bead
                        # apply to atom & neighbor
                        final[idx]=bead; final[n]=bead
                        merged=True; break
                if merged: break
            if merged: continue
            # try T->S convert
            for n in neigh:
                bead_str = final[n]
                if bead_str.startswith('T'):
                    new = 'S'+bead_str[1:]
                    ids=[i for i,v in enumerate(final) if v==bead_str]
                    for i in ids: final[i]=new
                    final[idx]=new; final[n]=new
                    merged=True; break
            if merged: continue
            # try S propagate
            for n in neigh:
                bead_str = final[n]
                if bead_str.startswith('S'):
                    ids=[i for i,v in enumerate(final) if v==bead_str]
                    for i in ids: final[i]=bead_str
                    final[idx]=bead_str; final[n]=bead_str
                    merged=True; break
            if merged: continue
            raise ValueError("Non-benzene 5-ring section not fully mappable!")
        elif final[idx]=='':
            raise ValueError("Non-benzene 5-ring section not fully mappable!")
    return final


# =============================================================================
# TASK 3.d: Determine if a Non-Ring Section Is 1-Bead Mappable
# =============================================================================
def is_non_ring_section_1bead_mappable(section: List[List[Any]]) -> bool:
    """
    Determine whether a non-ring section (section type == 0) is mappable as a single bead.
    
    Pseudocode:
      - Collect all edge atoms (atom[5] == True) into an array.
      - If exactly one edge atom exists, return True.
      - If none, throw an error.
      - If more than one, for each pair of edge atoms compute the shortest path (using BFS on the inner connections).
        If any such path is longer than 4 atoms, return False; otherwise, return True.
    """
    edge_atoms = [atom[0] for atom in section if atom[5]]
    if len(edge_atoms) == 1:
        return True
    if len(edge_atoms) == 0:
        raise ValueError("Cannot map non-ring section: no edge atoms found!")
    
    # Build a simple graph from inner connections
    graph = {}
    for atom in section:
        graph[atom[0]] = [section[t[0]][0] for t in atom[4]]
    
    def bfs_shortest_path(start: int, end: int) -> int:
        visited = set()
        queue = [(start, 1)]
        while queue:
            current, dist = queue.pop(0)
            if current == end:
                return dist
            visited.add(current)
            for neighbor in graph.get(current, []):
                if neighbor not in visited:
                    queue.append((neighbor, dist + 1))
        return float('inf')
    
    for i in range(len(edge_atoms)):
        for j in range(i+1, len(edge_atoms)):
            dist = bfs_shortest_path(edge_atoms[i], edge_atoms[j])
            if dist > 4:
                return False
    return True


# =============================================================================
# TASK 3.e: Map Non-Ring Section Using 1–Bead Mapping
# =============================================================================
def map_non_ring_section_1bead(section: List[List[Any]],
                               final: List[str],
                               martini_dict: Dict[str, List[Any]],
                               full_mapping: List[List[List[Any]]]) -> List[str]:
    """
    Map a non‐ring section (section type == 0) that is 1‐bead mappable.
    
    This implementation follows the pseudocode:
      - Build array0 as the list of global indices of edge atoms (atom[5] == True).
      - Depending on whether there are 2, 3, or 4 edge atoms, trace a SMILES‐like path
        (using inner connections from index 4) and then look for candidate bead keys in martini_dict.
      - For a 2-edge case, we trace a single linear path.
      - For a 3-edge case, we locate the atom with 3 inner connections and trace each branch.
      - For a 4-edge case, we trace each branch from the atom with 4 inner connections and then
        assign special bead names (“SX4e” to the branch edges and “TP1d” to the remaining atoms)
        if the unmapped atoms are a C–O pair.
        
    If exactly one candidate bead is found (or a special case is met) the bead name (with a random tag)
    is assigned to all atoms in the section. Otherwise, an error is raised.
    """
    # ----------------------------------------------------------------
    # Build cache: global index -> (section index, local index)
    global_to_sec_loc: Dict[int, (int,int)] = {}
    for sec_idx, sec in enumerate(full_mapping):
        for loc_idx, atom_info in enumerate(sec):
            gi = atom_info[0]
            global_to_sec_loc[gi] = (sec_idx, loc_idx)
    # ----------------------------------------------------------------
    # --- Helper functions ---
    def trace_linear_path(section: List[List[Any]], start_local: int) -> str:
        """Trace a linear SMILES path starting from an edge atom until a different edge is reached."""
        visited = set()
        current_local = start_local
        # Start with the atom’s type.
        path = section[current_local][1]
        prev_local = None
        visited.add(current_local)
        while True:
            # Get inner connections (list of tuples: (neighbor_local, bond order)).
            connections = section[current_local][4]
            next_local = None
            for (nbr, bond) in connections:
                if nbr != prev_local:
                    next_local = nbr
                    # Determine bond symbol.
                    bond_symbol = ""
                    if bond in (2, 1.5):
                        bond_symbol = "="
                    elif bond == 3:
                        bond_symbol = "#"
                    path += bond_symbol + section[nbr][1]
                    break
            if next_local is None:
                break
            # Stop if we reached a different edge atom.
            if section[next_local][5] and next_local != start_local:
                break
            if next_local in visited:
                break
            prev_local = current_local
            current_local = next_local
            visited.add(current_local)
        return path

    def trace_branch(section: List[List[Any]], center_local: int, next_local: int) -> str:
        """
        Trace a branch from a center atom along one inner connection.
        Returns a string representing that branch (using bond symbols and atom types).
        """
        visited = {center_local}
        # Find bond order from center to next_local.
        bond = None
        for (nbr, b) in section[center_local][4]:
            if nbr == next_local:
                bond = b
                break
        bond_symbol = ""
        if bond in (2, 1.5):
            bond_symbol = "="
        elif bond == 3:
            bond_symbol = "#"
        branch_str = bond_symbol + section[next_local][1]
        visited.add(next_local)
        current_local = next_local
        prev_local = center_local
        while True:
            connections = section[current_local][4]
            next_candidate = None
            for (nbr, b) in connections:
                if nbr != prev_local:
                    next_candidate = nbr
                    bond_sym = ""
                    if b in (2, 1.5):
                        bond_sym = "="
                    elif b == 3:
                        bond_sym = "#"
                    branch_str += bond_sym + section[nbr][1]
                    break
            if next_candidate is None:
                break
            if section[next_candidate][5]:
                break
            if next_candidate in visited:
                break
            prev_local = current_local
            current_local = next_candidate
            visited.add(current_local)
        return branch_str

    def get_final_edge(section: List[List[Any]], center_local: int, next_local: int) -> Optional[int]:
        """
        Trace along a branch (starting from center_local via next_local) and return the local index
        of the final edge atom (if found).
        """
        visited = {center_local}
        current_local = next_local
        prev_local = center_local
        visited.add(current_local)
        while True:
            connections = section[current_local][4]
            next_candidate = None
            for (nbr, _) in connections:
                if nbr != prev_local:
                    next_candidate = nbr
                    break
            if next_candidate is None:
                break
            if section[next_candidate][5]:
                return next_candidate
            if next_candidate in visited:
                break
            prev_local = current_local
            current_local = next_candidate
            visited.add(current_local)
        return current_local if section[current_local][5] else None

    # --- Main body of TASK 3.e ---
    # Build array0: list of global indices for edge atoms.
    edge_globals = [a[0] for a in section if a[5]]
    num_edges = len(edge_globals)
    # --- fallback for single-edge sections ---
    # … earlier in map_non_ring_section_1bead …
    # Build cache at top (you already have this)
    # global_to_sec_loc: Dict[int, (int,int)] = { … }
    if len(section) == 1:
        idx  = section[0][0]
        gi   = section[0][0]
        atom = section[0]

        # gather global neighbor indices (inner + outer)
        neigh = [ section[t[0]][0]           for t in atom[4] ] + \
                [ full_mapping[t[0]][t[1]][0] for t in atom[3] ]

        # Try merging into existing TC5/TC3 beads first
        for prefix in ("TC5","TC3"):
            for n in neigh:
                b = final[n]
                if b.startswith(prefix):
                    # clear the old group
                    old_tag = b
                    ids = [i for i,v in enumerate(final) if v == old_tag]
                    for i in ids:
                        final[i] = ""
                    # pick new S‐type key
                    bond = atom[4][0][1] if atom[4] else None
                    key  = pick_bead_key(martini_dict,
                                         atom[1].upper(),
                                         bond,
                                         kind='S')
                    if not key:
                        raise ValueError(f"Fallback: no S-key for {atom[1]} with bond {bond}")
                    bead = key + generate_random_string()
                    for i in ids:
                        final[i] = bead
                    final[gi] = final[n] = bead
                    return final

        # If it’s not C, and not O, we can’t handle it
        if atom[1].upper() != 'C':
            if atom[1].upper() != 'O':
                raise ValueError("1-edge non-ring not mappable (non-C)")
            
            # --- lone‐O special handling ---
            neighbor_beads = [final[n] for n in neigh]
            stripped = [
                (b[:-6].replace('+','') if b and len(b) > 6 else '')
                for b in neighbor_beads
            ]

            # Case 1: all neighbors unmapped
            if all(not s for s in stripped):
                raise ValueError("lone atom is next to unmapped atoms")

            # helper to get (sec,loc) for any global index
            def lookup_atom(gidx):
                sec_idx, loc_idx = global_to_sec_loc[gidx]
                return full_mapping[sec_idx][loc_idx]

            # Case 2: SC3
            if 'SC3' in stripped:
                i         = stripped.index('SC3')
                neigh_idx = neigh[i]
                bead_tag  = neighbor_beads[i]

                # try inner‐bond first
                bond = next((b for nbr, b in atom[4]
                             if section[nbr][0] == neigh_idx),
                            None)
                if bond is None:
                    # unpack the outer‐connection triples properly
                    bond = next((bo for sec_i, loc_i, bo in atom[3]
                                 if full_mapping[sec_i][loc_i][0] == neigh_idx),
                                None)

                # gather the existing SC3 group globals
                group_idxs = [j for j,v in enumerate(final) if v == bead_tag]

                # find the centre atom in that SC3 chain
                centre = None
                for g in group_idxs:
                    a = lookup_atom(g)
                    # count how many of its connections hit the other two SC3 atoms
                    # outer neighbours (already global):
                    outer_globals = [ full_mapping[si][li][0] for si,li,_ in a[3] ]
                    
                    # inner neighbours (lookup in a's own section):
                    sec_idx, _   = global_to_sec_loc[a[0]]
                    a_section    = full_mapping[sec_idx]
                    inner_globals = [ a_section[nbr_local][0] for nbr_local,_ in a[4] ]
                    
                    conns = outer_globals + inner_globals
                    if sum(1 for c in conns if c in group_idxs) == 2:
                        centre = a
                        break

                # is O connected to that centre?
                connected = False
                if centre:
                    connected = any(
                        section[t][0] == centre[0]
                        for t in atom[4]
                    )

                # pick new key
                if   bond == 1 and  connected: key = 'SN6'
                elif bond == 2 and  connected: key = 'SN6a'
                elif bond == 1 and not connected: key = 'N6'
                else:                            key = 'N6a'

                new_bead = key + generate_random_string()
                final[gi] = new_bead
                for g in group_idxs:
                    final[g] = new_bead
                return final

            # Case 3: SC4
            if 'SC4' in stripped:
                i         = stripped.index('SC4')
                neigh_idx = neigh[i]
                bead_tag  = neighbor_beads[i]

                # O–C bond
                bond = next((b for nbr,b in atom[4]
                             if section[nbr][0] == neigh_idx), None)
                if bond is None:
                    bond = next((bo for _,_,bo in atom[3]
                                 if lookup_atom(neigh_idx)[0] == neigh_idx),
                                None)

                key      = 'N6' if bond == 1 else 'N6a'
                new_bead = key + generate_random_string()

                final[gi] = new_bead
                # reassign the full SC4 group
                for j,v in enumerate(final):
                    if v == bead_tag:
                        final[j] = new_bead
                return final

            # Cases 4–6: TN4a, SN6, SN4a
            case_map = {
                'TN4a': {'prefixes': ['TN4a'], 'bond_map': {1: 'SN5a', 2: 'SN4a'}},
                'SN6':  {'prefixes': ['SN6'],  'bond_map': {1: 'N5',   2: 'N5a'}},
                'SN4a': {'prefixes': ['SN4a'], 'bond_map': {1: 'N5a',  2: 'N4a'}},
            }
            for base, info in case_map.items():
                if base in stripped:
                    i         = stripped.index(base)
                    bead_tag  = neighbor_beads[i]
                    neigh_idx = neigh[i]

                    ca = lookup_atom(neigh_idx)
                    if ca[1].upper() == 'C':
                        # find bond order: try inner then outer
                        bond = next((b for nbr,b in atom[4]
                                     if section[nbr][0] == neigh_idx),
                                    None)
                        if bond is None:
                            bond = next((bo for sec_i, loc_i, bo in atom[3]
                                         if full_mapping[sec_i][loc_i][0] == neigh_idx),
                                        None)

                        key = info['bond_map'].get(bond)
                        if key:
                            new_bead = key + generate_random_string()
                            # reassign the lone-O (gi) and the entire matched group
                            final[gi] = new_bead
                            for j, v in enumerate(final):
                                if v == bead_tag:
                                    final[j] = new_bead
                            return final
                    break  # if base matched but CA wasn’t C, skip to fallback

            # nothing matched
            print(f"[lone‑O] couldn’t map O@{gi}, neighbor beads = {neighbor_beads}")
            raise ValueError("1-edge non-ring not mappable (non-C)")
        # convert T* to S*
        for n in neigh:
            b=final[n]
            if b.startswith('T'):
                new='S'+b[1:]
                ids=[i for i,v in enumerate(final) if v==b]
                for i in ids: final[i]=new
                final[idx]=final[n]=new
                return final
        # propagate S*
        for n in neigh:
            b=final[n]
            if b.startswith('S'):
                ids=[i for i,v in enumerate(final) if v==b]
                for i in ids: final[i]=b
                final[idx]=final[n]=b
                return final
        raise ValueError("1-edge non-ring still not mappable")
      
    # CASE 1: 2 edge atoms.
    if num_edges == 2:
        start_global = edge_globals[0]
        start_local = next(i for i, atom in enumerate(section) if atom[0] == start_global)
        path_str = trace_linear_path(section, start_local).upper()

        # 1) Try section 7
        candidate_keys = [
            key for key, val in martini_dict.items()
            if val[0] == 7 and val[1] == 0
               and (val[2].upper() == path_str or val[2].upper() == path_str[::-1])
        ]
        candidate_keys = list(set(candidate_keys))

        # 2) Fallback to section 11 if needed
        if not candidate_keys:
            # ---- new: special split logic for 4‐atom sections ----
            # 2‑edge case
            if len(section) == 4:
                print(f"Linear path of 4 cannot be mapped: no candidate bead found for 2-edge mapping (path: {path_str}). Beginning to break down further.")
                # Build a little helper to map exactly two atoms by connectivity
                def map_two_atoms(local_a: int, local_b: int):
                    ga = section[local_a][0]
                    gb = section[local_b][0]
                    # Build the SMILES‐like path between them
                    sym = ""
                    for nbr, bo in section[local_a][4]:
                        if nbr == local_b:
                            sym = "=" if bo == 2 else "#" if bo == 3 else ""
                            break
                    path = (section[local_a][1] + sym + section[local_b][1]).upper()
            
                    # Lookup in section 7
                    keys = [
                        k for k,v in martini_dict.items()
                        if v[0] == 7 and v[1] == 0
                           and (v[2].upper() == path or v[2].upper() == path[::-1])
                    ]
                    if not keys:
                        # fallback to section 11
                        keys = [
                            k for k,v in martini_dict.items()
                            if v[0] == 11 and v[1] == 0
                               and (v[2].upper() == path or v[2].upper() == path[::-1])
                        ]
                    if not keys:
                        raise ValueError(
                            f"Linear path of 2 from 4 cannot be mapped: no candidate bead found for 2-edge mapping (path: {path})."
                        )
                    bead = keys[0] + generate_random_string()
                    final[ga] = bead
                    final[gb] = bead
            
                # Identify the first pair: start_global & its only neighbor
                start_local = next(i for i,a in enumerate(section) if a[0] == edge_globals[0])
                nbr_local, _ = section[start_local][4][0]
                # Identify the other two locals
                used = { section[start_local][0], section[nbr_local][0] }
                others = [i for i,a in enumerate(section) if a[0] not in used]
                if len(others) != 2:
                    raise ValueError("Expected exactly 2 leftover atoms when splitting 4‑atom 2‑edge case.")
                # Map each sub‑pair
                map_two_atoms(start_local, nbr_local)
                map_two_atoms(others[0],    others[1])
                return final

            else:
                candidate_keys = [
                    key for key, val in martini_dict.items()
                    if val[0] == 11 and val[1] == 0
                       and (val[2].upper() == path_str or val[2].upper() == path_str[::-1])
                ]
                candidate_keys = list(set(candidate_keys))
                if candidate_keys:
                    warnings.warn(
                        f"Warning: falling back to section 11 for 2-edge mapping (path: {path_str})",
                        UserWarning
                    )

        # 3) Commit or error
        if len(candidate_keys) == 1:
            rstr = generate_random_string()
            bead = candidate_keys[0] + rstr
            for atom in section:
                final[atom[0]] = bead
            return final

        # no candidates even after fallback
        if not candidate_keys:
            raise ValueError(
                f"Non-ring section cannot be mapped: no candidate bead found for 2-edge mapping (path: {path_str})."
            )
        else:
            # If more than one candidate remains, check for a special case.
            if path_str in ("CO", "OC"):
                for atom in section:
                    if atom[0] in edge_globals and atom[3]:
                        # Use the first outer connection, but instead of using final directly,
                        # look up the foreign section via full_mapping and take its first atom.
                        foreign_sec_index = atom[3][0][0]
                        index = atom[3][0][1]
                        foreign_atom = full_mapping[foreign_sec_index][index]  # the first atom in the foreign section
                        foreign_bead = final[foreign_atom[0]]
                        if foreign_bead.startswith("SX4e"):
                            rstr = generate_random_string()
                            for a in section:
                                final[a[0]] = candidate_keys[0] + rstr
                            return final
                        else:
                            rstr = generate_random_string()
                            for a in section:
                                final[a[0]] = candidate_keys[1] + rstr
                            return final
                    else:
                        rstr = generate_random_string()
                        for a in section:
                            final[a[0]] = candidate_keys[1] + rstr
                        return final
    # CASE 2: 3 edge atoms.
    elif num_edges == 3:
        # Find the atom with 3 inner connections.
        center_local = None
        for i, atom in enumerate(section):
            if len(atom[4]) == 3:
                center_local = i
                break
        if center_local is None:
            raise ValueError("Non‐ring section with 3 edges does not have an atom with 3 inner connections.")
        # Trace each branch from the center atom.
        branch_paths = []
        for (nbr, _) in section[center_local][4]:
            branch_str = trace_branch(section, center_local, nbr)
            branch_paths.append(branch_str)

        # 1) Try section 7 matching
        import re
        candidate_keys = []
        for key, val in martini_dict.items():
            if val[0] == 7 and val[1] == 1 and "(" in val[2]:
                prefix, remainder = val[2].split("(", 1)
                if prefix.strip() != section[center_local][1]:
                    continue
                branch_segments = re.findall(r'\(([^)]+)\)', val[2])
                if len(branch_segments) != len(branch_paths):
                    continue
                # order-independent match
                unmatched = branch_segments.copy()
                matched_all = True
                for bp in branch_paths:
                    for bs in list(unmatched):
                        if bp == bs:
                            unmatched.remove(bs)
                            break
                    else:
                        matched_all = False
                        break
                if matched_all:
                    candidate_keys.append(key)
        candidate_keys = list(set(candidate_keys))

        # 2) Fallback to section 11 if needed
        if not candidate_keys:
            for key, val in martini_dict.items():
                if val[0] == 11 and val[1] == 1 and "(" in val[2]:
                    prefix, remainder = val[2].split("(", 1)
                    if prefix.strip() != section[center_local][1]:
                        continue
                    branch_segments = re.findall(r'\(([^)]+)\)', val[2])
                    if len(branch_segments) != len(branch_paths):
                        continue
                    unmatched = branch_segments.copy()
                    matched_all = True
                    for bp in branch_paths:
                        for bs in list(unmatched):
                            if bp == bs:
                                unmatched.remove(bs)
                                break
                        else:
                            matched_all = False
                            break
                    if matched_all:
                        candidate_keys.append(key)
            if candidate_keys:
                warnings.warn(
                    f"Warning: falling back to section 11 for 3-edge mapping (center: {section[center_local][1]} and branches: {branch_paths})",
                    UserWarning
                )
            candidate_keys = list(set(candidate_keys))

        # 3) Commit or error
        if len(candidate_keys) == 1:
            rstr = generate_random_string()
            for atom in section:
                final[atom[0]] = candidate_keys[0] + rstr
            return final
        elif not candidate_keys:
            raise ValueError(
                f"Non‐ring section cannot be mapped: no candidate bead found for 3‐edge mapping (center: {section[center_local][1]} and branches: {branch_paths})."
            )
        else:
            raise ValueError(
                f"Non‐ring section cannot be mapped: ambiguous candidate keys for 3‐edge mapping (center: {section[center_local][1]} and branches: {branch_paths})."
            )
            ''' OLD VERSION KEPT JUST IN CASE
    # CASE 3: 4 edge atoms.
    elif num_edges == 4:
        # Find the center atom that has 4 inner connections.
        center_local = None
        for i, atom in enumerate(section):
            if len(atom[4]) == 4:
                center_local = i
                break
        if center_local is None:
            raise ValueError("Non‐ring section with 4 edges does not have an atom with 4 inner connections.")
        
        # Build array4: for each inner connection of the center atom, trace to the final edge atom.
        #  start by adding the central atom's global index.
        array4 = [section[center_local][0]]
        for (nbr, _) in section[center_local][4]:
            final_edge_local = get_final_edge(section, center_local, nbr)
            # Only add if the final atom is an edge and its type is 'F'
            if final_edge_local is not None and section[final_edge_local][5] and section[final_edge_local][1] == 'F':
                array4.append(section[final_edge_local][0])
        
        # Ensure exactly 4 branch edge atoms were found.
        if len(array4) != 4:
            raise ValueError("Non‐ring section cannot be mapped: expected 4 branch edge atoms, found {}.".format(len(array4)))
        
        # Assign the bead "SX4e" with a random tag to all atoms in array4.
        rstr = generate_random_string()
        for idx in array4:
            final[idx] = "SX4e" + rstr
        
        # Build array5: the list of atoms still unmapped.
        array5 = [atom[0] for atom in section if final[atom[0]] == ""]
        
        # If exactly 2 atoms remain and their types are C and O, assign bead "TP1d" with a random tag.
        if len(array5) == 2:
            types = [next(a for a in section if a[0] == idx)[1].upper() for idx in array5]
            if sorted(types) == ['C', 'O']:
                rstr2 = generate_random_string()
                for idx in array5:
                    final[idx] = "TP1d" + rstr2
                return final
            else:
                raise ValueError("Non‐ring section cannot be mapped: remaining atoms are not C and O (found: {}).".format(types))
        else:
            raise ValueError("Non‐ring section cannot be mapped: unexpected number of unmapped atoms in 4‐edge mapping (found {} in array5).".format(len(array5)))
    '''
    # CASE 3: 4 edge atoms.
    elif num_edges == 4:
        center_local = None
        for i, atom in enumerate(section):
            if len(atom[4]) == 4:
                center_local = i
                break
        if center_local is None:
            raise ValueError("Non‐ring section with 4 edges does not have an atom with 4 inner connections.")
        total_atoms = len(section)
        if total_atoms == 6:
            two_branch = None
            one_branches = []
            for (nbr_local, _) in section[center_local][4]:
                final_edge = get_final_edge(section, center_local, nbr_local)
                if final_edge is not None:
                    from collections import deque
        
                    def bfs_path(start, end):
                        visited = {start}
                        queue = deque([(start, [start])])
                        while queue:
                            curr, path = queue.popleft()
                            if curr == end:
                                return path
                            for (nbr, _) in section[curr][4]:
                                if nbr not in visited:
                                    visited.add(nbr)
                                    queue.append((nbr, path + [nbr]))
                        return []
        
                    path = bfs_path(center_local, final_edge)
                    if len(path) == 3:
                        two_branch = path   # [center_local, intermediate_local, final_edge]
                    else:
                        one_branches.append((center_local, final_edge))
        
            if two_branch is None or len(one_branches) != 3:
                raise ValueError("4-edge 6-atom split: could not identify branches correctly.")
        
            _, int_local, edge_local = two_branch
        
            # 1) Remove the inner-bond entry from center<->intermediate, but
            #    “move” it as an outer-connection on `int_local` pointing back to center.
            #    First, find the bond order between them:
            bond_to_center = None
            for (nbr, bo) in section[int_local][4]:
                if nbr == center_local:
                    bond_to_center = bo
                    break
            # Remove from int_local’s inner list:
            section[int_local][4] = [
                (nbr, bo) for (nbr, bo) in section[int_local][4]
                if nbr != center_local
            ]
            # Remove from center_local’s inner list:
            section[center_local][4] = [
                (nbr, bo) for (nbr, bo) in section[center_local][4]
                if nbr != int_local
            ]
        
            # Now “move” that edge into int_local’s outer list. We need the index of this section in full_mapping:
            sec_idx = full_mapping.index(section)
            section[int_local][3].append((sec_idx, center_local, bond_to_center))
        
            # 2) Build sub2 = [center_local] + three single leaves, recurse on sub2 first.
            sub2_locals = [center_local] + [loc for (_, loc) in one_branches]
            sub2 = [section[i] for i in sub2_locals]
            # Mark only the three leaf nodes as edges; center_local stays isedge=False.
            for (_, leaf_loc) in one_branches:
                section[leaf_loc][5] = True
            sub2_reindexed = reindex_section(sub2, sub2_locals)
            final = map_non_ring_section_1bead(sub2_reindexed, final, martini_dict, full_mapping)
        
            # 3) Now build sub1 = [int_local, edge_local], recurse on sub1.
            #    int_local already has [3] updated to list its outer = center.
            section[int_local][5] = True
            section[edge_local][5] = True
            sub1_locals = [int_local, edge_local]
            sub1 = [section[i] for i in sub1_locals]
            sub1_reindexed = reindex_section(sub1, sub1_locals)
            final = map_non_ring_section_1bead(sub1_reindexed, final, martini_dict, full_mapping)
        
            return final

        # Part 2: if section has 5 atoms, treat like 3-edge but with val[1]==2
        elif total_atoms == 5:
            # All four neighbors of center are direct edges
            branch_paths = []
            for (nbr, _) in section[center_local][4]:
                if section[nbr][5]:
                    branch_str = trace_branch(section, center_local, nbr)
                    branch_paths.append(branch_str)
            if len(branch_paths) != 4:
                raise ValueError(f"Non‐ring 4-edge 5-atom mapping: expected 4 branches, found {len(branch_paths)}.")
            import re
            candidate_keys = []
            for key, val in martini_dict.items():
                if val[0] == 7 and val[1] == 2 and "(" in val[2]:
                    prefix, remainder = val[2].split("(", 1)
                    if prefix.strip() != section[center_local][1]:
                        continue
                    branch_segments = re.findall(r'\(([^)]+)\)', val[2])
                    if len(branch_segments) != len(branch_paths):
                        continue
                    unmatched = branch_segments.copy()
                    matched_all = True
                    for bp in branch_paths:
                        for bs in list(unmatched):
                            if bp == bs:
                                unmatched.remove(bs)
                                break
                        else:
                            matched_all = False
                            continue
                    if val[4] == "." and len(unmatched) <= 2:
                        candidate_keys.append(key)
                    if matched_all:
                        candidate_keys.append(key)
            candidate_keys = list(set(candidate_keys))
            if len(candidate_keys) == 1:
                rstr = generate_random_string()
                for atom in section:
                    final[atom[0]] = candidate_keys[0] + rstr
                return final
            elif not candidate_keys:
                raise ValueError(
                    f"Non‐ring section cannot be mapped: no candidate bead found for 4-edge mapping (center: {section[center_local][1]} and branches: {branch_paths})."
                )
            else:
                raise ValueError(
                    f"Non‐ring section cannot be mapped: ambiguous candidate keys for 4-edge mapping (center: {section[center_local][1]} and branches: {branch_paths})."
                )
        else:
            raise ValueError(f"Non‐ring 4-edge mapping only supported for 5 or 6 atoms; found {total_atoms}.")

    # CASE 4: 5 or 6 edge atoms.
    elif num_edges > 4:
        # find the two centers (atoms with >2 inner connections)
        centers = [i for i, atom in enumerate(section) if len(atom[4]) > 2]
        if len(centers) != 2:
            raise ValueError(f"Need exactly 2 centers to split; found {len(centers)}")
    
        c0, c1 = centers
    
        # disconnect the two centers in this section
        section[c0][4] = [(nbr, bo) for nbr, bo in section[c0][4] if nbr != c1]
        section[c1][4] = [(nbr, bo) for nbr, bo in section[c1][4] if nbr != c0]
    
        # make a local copy of the mapping so we don't clobber the caller's
        local_mapping = full_mapping[:]
        try:
            sec_idx = local_mapping.index(section)
        except ValueError:
            sec_idx = None
    
        # trace out each sub‑section by DFS (they’re now disconnected)
        sub_locals_list = []
        for center in centers:
            visited = {center}
            stack   = [center]
            sub_locs = []
            while stack:
                cur = stack.pop()
                sub_locs.append(cur)
                for nbr, _ in section[cur][4]:
                    if nbr not in visited:
                        visited.add(nbr)
                        stack.append(nbr)
            sub_locals_list.append(sorted(sub_locs))
    
        # rebuild each sub‑section with correct local indices
        subs = []
        for locs in sub_locals_list:
            # extract the raw atoms
            raw = [section[i] for i in locs]
            # reindex so raw[*][4] refers to positions within this small list
            sub = reindex_section(raw, locs)
            # recompute isedge flags
            fix_isedge(sub)
            subs.append(sub)
    
        # splice them into our local mapping (if we know where this section lived)
        if sec_idx is not None:
            local_mapping.pop(sec_idx)
            # insert in reverse so subs[0] ends up at sec_idx
            local_mapping.insert(sec_idx, subs[1])
            local_mapping.insert(sec_idx, subs[0])
    
        # recurse on each piece using our updated local_mapping
        for sub_section in subs:
            final = map_non_ring_section_1bead(
                sub_section,
                final,
                martini_dict,
                local_mapping
            )
        return final
    # Should not be reached.
    print(final)
    raise ValueError("Non‐ring section cannot be mapped: unhandled edge count case.")
    return final


# =============================================================================
# TASK 3.f: Map Non-Ring Sections That Are Too Long (Subdivision)
# =============================================================================
def fix_isedge(section: List[List[Any]]) -> None:
    """
    Recalculate the isedge flag for each atom in a section.
    An atom is considered an edge if it has 0 or 1 inner connections (index 4) within the section.
    """
    for atom in section:
        if len(atom[4]) <= 1:
            atom[5] = True
        else:
            atom[5] = False


def subdivide_non_ring_section_normal(section: List[List[Any]], 
                                      final: List[str],
                                      martini_dict: Dict[str, List[Any]],
                                      full_mapping: List[List[List[Any]]]) -> Tuple[List[List[Any]], List[List[Any]], List[int], List[int]]:
    """
    Normal subdivision using DFS from a single candidate edge.
    Returns:
      sub_section, remainder, sub_old_indices, rem_old_indices
    """
    # Build candidate_edges using local indices (atoms flagged as edge).
    candidate_edges = [i for i, atom in enumerate(section) if atom[5]]
    if not candidate_edges:
        candidate_edges = [0, len(section) - 1]
        
    # SPECIAL CASE for linear sections (exactly 2 edge nodes)
    if len(candidate_edges) == 2:
        if len(section) == 9:
            sub_old_indices = list(range(3))
            rem_old_indices = list(range(3, 9))
            sub_section = [section[i] for i in sub_old_indices]
            remainder = [section[i] for i in rem_old_indices]
            return sub_section, remainder, sub_old_indices, rem_old_indices
        elif len(section) == 6:
            sub_old_indices = list(range(3))
            rem_old_indices = list(range(3, 6))
            sub_section = [section[i] for i in sub_old_indices]
            remainder = [section[i] for i in rem_old_indices]
            return sub_section, remainder, sub_old_indices, rem_old_indices
        
    # (Inside subdivide_non_ring_section_normal)
    if len(candidate_edges) > 2:
        from collections import deque

        # Build a graph of inner connections (using local indices).
        graph = {i: [tup[0] for tup in atom[4]] for i, atom in enumerate(section)}

        def bfs_local(start: int, end: int) -> int:
            visited = set()
            queue = deque([(start, 0)])
            while queue:
                current, dist = queue.popleft()
                if current == end:
                    return dist
                visited.add(current)
                for nbr in graph.get(current, []):
                    if nbr not in visited:
                        queue.append((nbr, dist + 1))
            return float('inf')

        # Make a copy of candidate_edges before removal.
        candidate_edges_copy = candidate_edges[:]

        # Remove too-close candidate edges.
        to_remove = set()
        for i in range(len(candidate_edges)):
            for j in range(i + 1, len(candidate_edges)):
                if bfs_local(candidate_edges[i], candidate_edges[j]) <= 4:
                    to_remove.add(candidate_edges[i])
                    to_remove.add(candidate_edges[j])
        candidate_edges = [i for i in candidate_edges if i not in to_remove]

        if not candidate_edges:
            # --- New procedure starts here ---
            # Reset candidate_edges from the copy.
            candidate_edges = candidate_edges_copy

            # Build an array of candidate edge pairs with their BFS distances.
            minimal_pairs = []
            min_distance = float('inf')
            for i in range(len(candidate_edges)):
                for j in range(i + 1, len(candidate_edges)):
                    d = bfs_local(candidate_edges[i], candidate_edges[j])
                    if d < min_distance:
                        min_distance = d
                        minimal_pairs = [[candidate_edges[i], d, candidate_edges[j]]]
                    elif d == min_distance:
                        pair_sorted = sorted([candidate_edges[i], candidate_edges[j]])
                        # Avoid duplicates (regardless of order)
                        if pair_sorted not in [sorted([p[0], p[2]]) for p in minimal_pairs]:
                            minimal_pairs.append([candidate_edges[i], d, candidate_edges[j]])

            def bfs_path_between(start: int, end: int, section: List[List[Any]]) -> List[int]:
                """Return a BFS path from start to end as a list of local indices, or empty if none."""
                q = deque([(start, [start])])
                visited = {start}
                while q:
                    current, path = q.popleft()
                    if current == end:
                        return path
                    for (nbr, bond) in section[current][4]:
                        if nbr not in visited:
                            visited.add(nbr)
                            q.append((nbr, path + [nbr]))
                return []

            def is_connected(indices: List[int], section: List[List[Any]]) -> bool:
                """Check if the atoms with local indices in 'indices' form a connected subgraph."""
                if not indices:
                    return False
                visited = set()
                to_visit = [indices[0]]
                while to_visit:
                    curr = to_visit.pop(0)
                    visited.add(curr)
                    for (nbr, _) in section[curr][4]:
                        if nbr in indices and nbr not in visited:
                            to_visit.append(nbr)
                return set(indices) == visited

            def get_connected_components(indices: List[int], section: List[List[Any]]) -> List[List[int]]:
                """Return connected components (lists of indices) in the subgraph induced by 'indices'."""
                components = []
                remaining = set(indices)
                while remaining:
                    start = next(iter(remaining))
                    comp = []
                    q = [start]
                    while q:
                        curr = q.pop(0)
                        if curr not in comp:
                            comp.append(curr)
                            for (nbr, _) in section[curr][4]:
                                if nbr in remaining and nbr not in comp:
                                    q.append(nbr)
                    components.append(comp)
                    remaining -= set(comp)
                return components

            # Now try each candidate minimal pair.
            success = False
            chosen_sub_old_indices = None
            chosen_rem_old_indices = None
            for pair_info in minimal_pairs:
                edge1, d, edge2 = pair_info
                path = bfs_path_between(edge1, edge2, section)
                if not path:
                    continue
                candidate_sub_old_indices = sorted(path)
                candidate_rem_old_indices = [i for i in range(len(section)) if i not in candidate_sub_old_indices]
                # First, ensure connectivity on the remainder:
                if not is_connected(candidate_rem_old_indices, section):
                    section_index = full_mapping.index(section)
                    # The remainder is split into multiple connected components.
                    components = get_connected_components(candidate_rem_old_indices, section)
                    if len(components) == 2:
                        if len(components[0]) == len(components[1]):
                            # Try both merging options.
                            option1_sub = sorted(candidate_sub_old_indices + components[0])
                            option1_rem = [i for i in range(len(section)) if i not in option1_sub]
                            option2_sub = sorted(candidate_sub_old_indices + components[1])
                            option2_rem = [i for i in range(len(section)) if i not in option2_sub]
                            try:
                                dummy_final = final[:]  # Test option1
                                # Automatically find the section index in full_mapping.
                                section_index = full_mapping.index(section)
                                dummy_full_mapping = full_mapping[:]  
                                # Remove the section being subdivided
                                dummy_full_mapping.pop(section_index)
                                # Insert the candidate subdivisions in its place.
                                dummy_full_mapping.insert(section_index, [section[i] for i in option1_sub])
                                dummy_full_mapping.insert(section_index + 1, [section[i] for i in option1_rem])
                                # --- New reindexing step:
                                reindexed_sub = reindex_section([section[i] for i in option1_sub], option1_sub)
                                reindexed_rem = reindex_section([section[i] for i in option1_rem], option1_rem)
                                # Now update the boundary mapping on the reindexed candidate subdivisions.
                                updated_sub_section, updated_remainder = update_divided_section_mapping(
                                    reindexed_sub,
                                    reindexed_rem
                                )
                                dummy_final = map_non_ring_section_1bead(updated_sub_section, dummy_final, martini_dict, dummy_full_mapping)
                                dummy_final = map_non_ring_section_1bead(updated_remainder, dummy_final, martini_dict, dummy_full_mapping)
                                candidate_sub_old_indices = option1_sub
                                candidate_rem_old_indices = option1_rem
                            except Exception:
                                try:
                                    dummy_final = final[:]  # Test option2
                                    section_index = full_mapping.index(section)
                                    dummy_full_mapping = full_mapping[:]
                                    dummy_full_mapping.pop(section_index)
                                    dummy_full_mapping.insert(section_index, [section[i] for i in option2_sub])
                                    dummy_full_mapping.insert(section_index + 1, [section[i] for i in option2_rem])
                                    # --- New reindexing step:
                                    reindexed_sub = reindex_section([section[i] for i in option2_sub], option2_sub)
                                    reindexed_rem = reindex_section([section[i] for i in option2_rem], option2_rem)
                                    updated_sub_section, updated_remainder = update_divided_section_mapping(
                                        reindexed_sub,
                                        reindexed_rem
                                    )
                                    dummy_final = map_non_ring_section_1bead(updated_sub_section, dummy_final, martini_dict, dummy_full_mapping)
                                    dummy_final = map_non_ring_section_1bead(updated_remainder, dummy_final, martini_dict, dummy_full_mapping)
                                    candidate_sub_old_indices = option2_sub
                                    candidate_rem_old_indices = option2_rem
                                except Exception:
                                    continue
                        else:
                            # When the two components have different sizes, choose the smaller one.
                            if len(components[0]) > len(components[1]):
                                merged_sub = sorted(candidate_sub_old_indices + components[1])
                            else:
                                merged_sub = sorted(candidate_sub_old_indices + components[0])
                            rem_after_merge = [i for i in range(len(section)) if i not in merged_sub]
                            if not is_connected(rem_after_merge, section):
                                continue
                            candidate_sub_old_indices = merged_sub
                            candidate_rem_old_indices = rem_after_merge
                    else:
                        continue
                # *** New test: update boundary mapping and try mapping the candidate parts.
                try:
                    dummy_final = final[:]  # Make a copy of current final mapping.
                    section_index = full_mapping.index(section)
                    dummy_full_mapping = full_mapping[:]  # Copy the full mapping.
                    dummy_full_mapping.pop(section_index)
                    dummy_full_mapping.insert(section_index, [section[i] for i in candidate_sub_old_indices])
                    dummy_full_mapping.insert(section_index + 1, [section[i] for i in candidate_rem_old_indices])
                    # --- New reindexing step:
                    reindexed_sub = reindex_section([section[i] for i in candidate_sub_old_indices], candidate_sub_old_indices)
                    reindexed_rem = reindex_section([section[i] for i in candidate_rem_old_indices], candidate_rem_old_indices)
                    updated_sub_section, updated_remainder = update_divided_section_mapping(
                        reindexed_sub,
                        reindexed_rem
                    )
                    dummy_final = map_non_ring_section_1bead(updated_sub_section, dummy_final, martini_dict, dummy_full_mapping)
                    dummy_final = map_non_ring_section_1bead(updated_remainder, dummy_final, martini_dict, dummy_full_mapping)
                except Exception:
                    # Mapping candidate failed; try next minimal pair.
                    continue

                # If we reached this point without error, candidate mapping succeeded.
                chosen_sub_old_indices = candidate_sub_old_indices
                chosen_rem_old_indices = candidate_rem_old_indices
                success = True
                break
            if success and chosen_sub_old_indices is not None and chosen_rem_old_indices is not None:
                # Instead of proceeding with DFS tracing, immediately return this subdivision.
                sub_section = [section[i] for i in chosen_sub_old_indices]
                remainder = [section[i] for i in chosen_rem_old_indices]
                return sub_section, remainder, chosen_sub_old_indices, chosen_rem_old_indices
            else:
                raise ValueError("Non‐ring section division leads to unassignable sections.")
        # --- End of new procedure ---
    # Use the first candidate as seed for DFS tracing if candidate_edges is nonempty.
    seed = candidate_edges[0]

    def dfs_trace(local_idx: int, visited: set) -> List[int]:
        trace = [local_idx]
        visited.add(local_idx)
        for (nbr, bond) in section[local_idx][4]:
            if nbr not in visited and len(section[nbr][4]) < 3:
                trace.extend(dfs_trace(nbr, visited))
                break
        return trace

    visited_trace = set()
    trace_nodes = dfs_trace(seed, visited_trace)
    # NEW CONDITION 1: If the trace is long enough and the fourth node has ≥ 3 inner connections, limit trace to 3 nodes.
    if len(trace_nodes) >= 4 and len(section[trace_nodes[3]][4]) >= 3:
        trace_nodes = trace_nodes[:3]
    if len(trace_nodes) > 4:
        trace_nodes = trace_nodes[:4]
    # NEW CONDITION 2: Ensure remainder has at least 2 atoms.
    while len(section) - len(trace_nodes) < 2 and len(trace_nodes) > 2:
        trace_nodes = trace_nodes[:-1]
    if len(section) - len(trace_nodes) < 2:
        raise ValueError("Subdivision not possible: remainder too small (normal method).")
    sub_old_indices = sorted(trace_nodes)
    rem_old_indices = [i for i in range(len(section)) if i not in sub_old_indices]
    sub_section = [section[i] for i in sub_old_indices]
    remainder = [section[i] for i in rem_old_indices]
    return sub_section, remainder, sub_old_indices, rem_old_indices


def subdivide_non_ring_section_multi(section: List[List[Any]]) -> Tuple[List[List[Any]], List[List[Any]], List[int], List[int]]:
    """
    Special subdivision: when candidate_edges has 2 or more elements.
    Find the shortest path between any two edge nodes (without the <=4 constraint),
    and use that path as the sub_section.
    Returns:
      sub_section, remainder, sub_old_indices, rem_old_indices
    """
    from collections import deque
    edge_nodes = [i for i, atom in enumerate(section) if atom[5]]
    if len(edge_nodes) < 2:
        raise ValueError("Special subdivision not possible: fewer than 2 edge nodes.")
    graph = {i: [tup[0] for tup in atom[4]] for i, atom in enumerate(section)}
    def bfs_path(start: int, end: int) -> Tuple[int, List[int]]:
        queue = deque([(start, [start])])
        visited = {start}
        while queue:
            current, path = queue.popleft()
            if current == end:
                return (len(path)-1, path)
            for nbr in graph.get(current, []):
                if nbr not in visited:
                    visited.add(nbr)
                    queue.append((nbr, path + [nbr]))
        return (float('inf'), [])
    best_path = None
    best_length = float('inf')
    for i in range(len(edge_nodes)):
        for j in range(i+1, len(edge_nodes)):
            length, path = bfs_path(edge_nodes[i], edge_nodes[j])
            if length < best_length:
                best_length = length
                best_path = path
    if best_path is None or len(best_path) < 2:
        raise ValueError("Special subdivision not possible: no valid path found.")
    sub_old_indices = sorted(best_path)
    rem_old_indices = [i for i in range(len(section)) if i not in sub_old_indices]
    if len(rem_old_indices) < 2:
        # If remainder too small, try removing the last node from best_path.
        sub_old_indices = best_path[:-1]
        rem_old_indices = [i for i in range(len(section)) if i not in sub_old_indices]
        if len(rem_old_indices) < 2:
            raise ValueError("Special subdivision not possible: remainder too small.")
    sub_section = [section[i] for i in sub_old_indices]
    remainder = [section[i] for i in rem_old_indices]
    return sub_section, remainder, sub_old_indices, rem_old_indices


def subdivide_non_ring_section(section: List[List[Any]],
                               final: List[str],
                               full_mapping: List[List[List[Any]]],
                               martini_dict: Dict[str, List[Any]]) -> Tuple[List[List[Any]], List[List[Any]], List[int], List[int]]:
    """
    Chooses the subdivision method based on candidate edges.
    If candidate_edges has <=1 element, use normal; if >=2, try normal first and fall back to multi-edge.
    """
    candidate_edges = [i for i, atom in enumerate(section) if atom[5]]
    if len(candidate_edges) <= 1:
        return subdivide_non_ring_section_normal(section, final, martini_dict, full_mapping)
    else:
        try:
            return subdivide_non_ring_section_normal(section, final, martini_dict, full_mapping)
        except ValueError as e:
            print("Critical error:", e)
            return subdivide_non_ring_section_multi(section)


def reindex_section(section: List[List[Any]], old_indices: List[int]) -> List[List[Any]]:
    """
    Reindex a subdivided section so that inner connections (index 4) refer to new local indices.
    'old_indices' is a list of the original positions (local indices) for the atoms in this new section.
    """
    mapping = {old: new for new, old in enumerate(old_indices)}
    new_section = []
    for new_idx, atom in enumerate(section):
        new_atom = list(atom)
        new_inner = []
        for (nbr_old, bond) in new_atom[4]:
            # Convert neighbor's old index to new local index if it appears in the mapping.
            if nbr_old in mapping:
                new_inner.append((mapping[nbr_old], bond))
        new_atom[4] = new_inner
        new_section.append(new_atom)
    return new_section


def update_divided_section_mapping(sub_section: List[List[Any]], remainder: List[List[Any]]
                                  ) -> Tuple[List[List[Any]], List[List[Any]]]:
    """
    TASK 3.g (Normal Version, simplified):
    Given a subdivision into two parts (sub_section and remainder) along with their original positions
    (sub_old_indices and rem_old_indices, respectively),
    update the mapping so that the inner connections are reindexed to refer to the new local positions.
    Then update the isedge flag for both parts.
    
    Outer connections (index 3) are left unchanged.
    """
    fix_isedge(sub_section)
    fix_isedge(remainder)
    return sub_section, remainder

def map_non_ring_section_long(section: List[List[Any]],
                              final: List[str],
                              martini_dict: Dict[str, List[Any]],
                              full_mapping: List[List[List[Any]]],
                              section_index: int) -> Tuple[List[str], List[List[List[Any]]]]:
    """
    TASK 3.f/3.g (Revised Special):
    For a non-ring section (type == 0) that is too long to be mapped as a single bead,
    subdivide the section (using either normal or multi-edge splitting as needed) and update the mapping.
    
    Process:
      - Subdivide the section using subdivide_non_ring_section.
      - Reindex each new section.
      - Depending on which method was used, update the boundary.
        (For multi-edge, use update_divided_section_mapping_multi; otherwise, use normal.)
      - Replace the section in full_mapping with the two new sections.
      - Recursively map each new section.
    Returns updated final and full mapping.
    """
    sub_sec, rem_sec, sub_old_indices, rem_old_indices = subdivide_non_ring_section(section, final, full_mapping, martini_dict)
    sub_sec = reindex_section(sub_sec, sub_old_indices)
    rem_sec = reindex_section(rem_sec, rem_old_indices)
    
    sub_sec, rem_sec = update_divided_section_mapping(sub_sec, rem_sec)
    
    full_mapping.pop(section_index)
    full_mapping.insert(section_index, sub_sec)
    full_mapping.insert(section_index + 1, rem_sec)
    
    final = map_non_ring_section_1bead(sub_sec, final, martini_dict, full_mapping)
    if not is_non_ring_section_1bead_mappable(rem_sec):
        final, full_mapping = map_non_ring_section_long(rem_sec, final, martini_dict, full_mapping, section_index + 1)
    else:
        final = map_non_ring_section_1bead(rem_sec, final, martini_dict, full_mapping)
    return final, full_mapping

# =============================================================================
# TASK 3.h.a: Map Non-Benzene 3-Ring Section
# =============================================================================
def map_nonbenzene_3_ring_section(section: List[List[Any]],
                                  final: List[str],
                                  martini_dict: Dict[str, List[Any]]) -> List[str]:
    """
    TASK 3.h.a: Map a non-benzene 3-membered ring section.
    
    For now, trace the atoms in the section and if they are all carbon atoms 
    with exactly two inner connections and every inner bond has order 1,
    assign the bead "SC3" + random string to each atom in the section.
    Otherwise, throw an error indicating that the 3-ring is not mappable.
    """
    # Check that every atom is carbon, has exactly 2 inner connections, and that each bond is 1.
    for atom in section:
        # Verify atom is carbon.
        if atom[1].upper() != "C":
            raise ValueError("3-ring not mappable: Found an atom that is not carbon.")
        # Verify there are exactly 2 inner connections.
        if len(atom[4]) != 2:
            raise ValueError("3-ring not mappable: Expected exactly 2 inner connections per atom.")
        # Verify that each inner connection has a bond order of 1.
        for (nbr, bond) in atom[4]:
            if bond != 1:
                raise ValueError("3-ring not mappable: Found an inner bond order not equal to 1.")
    
    # All atoms passed the checks.
    bead = "SC3" + generate_random_string()
    for atom in section:
        final[atom[0]] = bead
    return final

# =============================================================================
# TASK 3.h.b: Map Non-Benzene 4-Ring Section
# =============================================================================
def map_nonbenzene_4_ring_section(section: List[List[Any]],
                                  final: List[str],
                                  martini_dict: Dict[str, List[Any]]) -> List[str]:
    raise ValueError("Non-benzene ring section of size 4 is not supported.")


# =============================================================================
# Top-Level Mapping Function (Tree-Like Algorithm)
# =============================================================================
def map_martini_beads(mapping: List[List[List[Any]]],
                      final: List[str],
                      martini_dict: Dict[str, List[Any]]) -> List[str]:
    """
    Two-pass ring handling + then non-ring sections.
    """
    skipped = []

    # --- PASS 1: Benzene rings (type 2) ---
    for idx, section in enumerate(mapping):
        if section and section[0][2] == 2 and any(final[a[0]] == "" for a in section):
            try:
                new_final = map_benzene_ring_section(section, final.copy(), martini_dict, mapping)
                final = new_final
            except Exception as e:
                # record for second pass
                print(f"mapping section {idx} failed:", e)
                skipped.append(("benzene", idx))

    # --- PASS 1: Non-benzene rings (type 1) ---
    for idx, section in enumerate(mapping):
        if section and section[0][2] == 1 and any(final[a[0]] == "" for a in section):
            size = len(section)
            try:
                if size == 6:
                    new_final = map_nonbenzene_6_ring_section(section, final.copy(), martini_dict, mapping)
                elif size == 5:
                    try:
                        new_final = map_nonbenzene_5_ring_section(section, final.copy(), martini_dict, mapping)
                    except Exception as e:
                        print(f"mapping section {idx} failed as a 5-ring section, attempt as if it is a 6-ring section:", e)
                        new_final = map_nonbenzene_6_ring_section(section, final.copy(), martini_dict, mapping)
                elif size == 4:
                    new_final = map_nonbenzene_6_ring_section(section, final.copy(), martini_dict, mapping)
                elif size == 3:
                    new_final = map_nonbenzene_3_ring_section(section, final.copy(), martini_dict)
                else:
                    raise ValueError(f"Unexpected non-benzene ring size {size}")
                final = new_final
            except Exception as e:
                print(f"mapping section {idx} failed:", e)
                skipped.append(("nonbenzene", idx))

    # --- PASS 2: Retry skipped ---
    for sect_type, idx in skipped:
        section = mapping[idx]
        if sect_type == "benzene":
            final = map_benzene_ring_section(section, final, martini_dict, mapping)
        else:
            size = len(section)
            if size == 6:
                final = map_nonbenzene_6_ring_section(section, final, martini_dict, mapping)
            elif size == 5:
                try:
                    new_final = map_nonbenzene_5_ring_section(section, final.copy(), martini_dict, mapping)
                    final = new_final
                except Exception as e:
                    print(f"mapping section {idx} failed as a 5-ring section, attempt as if it is a 6-ring section:", e)
                    new_final = map_nonbenzene_6_ring_section(section, final.copy(), martini_dict, mapping)
                    final = new_final
            elif size == 4:
                new_final = map_nonbenzene_6_ring_section(section, final.copy(), martini_dict, mapping)
            elif size == 3:
                final = map_nonbenzene_3_ring_section(section, final, martini_dict)
    # --- Then the non-ring sections as before ---
    i = 0
    while i < len(mapping):
        section = mapping[i]
        if section and section[0][2] == 0:
            # keep going until every atom in this section is filled
            while any(final[a[0]] == "" for a in section):
                if is_non_ring_section_1bead_mappable(section):
                    final = map_non_ring_section_1bead(section, final, martini_dict, mapping)
                else:
                    final, mapping = map_non_ring_section_long(section, final, martini_dict, mapping, i)
                    section = mapping[i]
        i += 1
    # sanity check
    if any(name == "" for name in final):
        raise ValueError("Final mapping incomplete; some atoms remain unmapped!")
    return final
