from output_manip import get_coordinates_from_smiles
import math
def group_beads_by_type(final):
    """
    Group indices from the final list by unique bead string.
    For each unique bead string (in order of first appearance), all indices where
    that string occurs are collected into one group.
    
    Returns:
      - groups: list of lists, where each sublist contains the indices (from final)
                that belong to one unique bead.
      - bead_types: list of the corresponding bead strings (unfixed) for each group.
    """
    groups_dict = {}
    order = []
    for i, bead in enumerate(final):
        if bead not in groups_dict:
            groups_dict[bead] = []
            order.append(bead)
        groups_dict[bead].append(i)
    groups = []
    bead_types = []
    for bead in order:
        groups.append(groups_dict[bead])
        bead_types.append(bead)
    return groups, bead_types

def find_atom_in_mapping(j, mapping):
    """
    Search through every section in mapping to find the first atom where atom[0] == j.
    Returns a tuple (section_index, section, atom) or (None, None, None) if not found.
    """
    for sec_index, section in enumerate(mapping):
        for atom in section:
            if atom[0] == j:
                return sec_index, section, atom
    return None, None, None

def get_composed_elements(groups, mapping):
    """
    For each bead group (list of indices), search through every section in mapping
    for the atom whose trace index (atom[0]) matches and append the element (atom[1]).
    
    Returns:
      - coe: list of lists where each inner list contains the element for every atom in that bead.
    """
    coe = []
    for group in groups:
        elements = []
        for j in group:
            _, _, atom = find_atom_in_mapping(j, mapping)
            if atom is not None:
                elements.append(atom[1])
        coe.append(elements)
    return coe

def compute_beads_connections(coi, mapping):
    """
    For each bead (given by coi) and for every atom in that bead, search every section in mapping
    for the atom with trace index j. Then:
      - For each neighbor in the internal connections (atom[4]), the tuple is 
        (local_index, bond). The actual trace index is obtained from the same section:
        section[local_index][0].
      - For each connection in the external connections (atom[3]), the tuple is 
        [section_index, local_index, bond]. The actual trace index is obtained as
        mapping[section_index][local_index][0].
        
    If the retrieved trace index is not in the current bead, then find which bead group contains it
    and record the connection.
    
    To avoid duplicates, each connection is converted to a canonical form (i.e. (i, j, bond) with i <= j)
    and then sorted.
    
    Returns three lists: (origin_beads, connected_beads, bond_types)
    """
    connections = []
    
    # For each bead group i.
    for i, group in enumerate(coi):
        for j in group:
            sec_index, section, atom = find_atom_in_mapping(j, mapping)
            if atom is None:
                continue
            # Check internal connections (index 4)
            for neighbor in atom[4]:
                local_index = neighbor[0]   # location in the same section
                bond = neighbor[1]
                # Retrieve the actual trace index from the current section.
                neighbor_trace = section[local_index][0]
                if neighbor_trace not in group:
                    for c_index, other_group in enumerate(coi):
                        if neighbor_trace in other_group:
                            connections.append((i, c_index, bond))
                            break
            # Check external connections (index 3)
            for ext in atom[3]:
                # ext: [section_index, local_index, bond]
                foreign_sec_index = ext[0]
                local_index = ext[1]
                bond = ext[2]
                foreign_trace = mapping[foreign_sec_index][local_index][0]
                if foreign_trace not in group:
                    for c_index, other_group in enumerate(coi):
                        if foreign_trace in other_group:
                            connections.append((i, c_index, bond))
                            break
    # Remove duplicates by enforcing canonical order (i <= j).
    unique_connections = set()
    for (i, j, bond) in connections:
        if i > j:
            i, j = j, i
        unique_connections.add((i, j, bond))
    # Sort the unique connections by origin then connection index.
    sorted_connections = sorted(unique_connections, key=lambda x: (x[0], x[1]))
    origin_beads = [conn[0] for conn in sorted_connections]
    connected_beads = [conn[1] for conn in sorted_connections]
    bond_types = [conn[2] for conn in sorted_connections]
    return origin_beads, connected_beads, bond_types

def compute_num_connections(num_beads, origin_beads, connected_beads):
    """
    Count the number of connections for each bead symmetrically.
    """
    num_connections = [0] * num_beads
    for bead_origin, bead_connected in zip(origin_beads, connected_beads):
        num_connections[bead_origin] += 1
        num_connections[bead_connected] += 1  # increment both beads symmetrically
    return num_connections

def fix_beadtypes(bead_types):
    """
    For every bead type string, remove the last 6 characters and any '+' symbols.
    (Adjust as needed based on the actual input format.)
    """
    fixed = []
    for b in bead_types:
        fixed_b = b[:-6].replace('+', '')
        fixed.append(fixed_b)
    return fixed

def export_bead_mapping(final, mapping, smiles, compound_name, write_file=True):
    """
    Build bead‐groupings and connections from your Martini‐3 mapping,
    compute centroids & masses, and write a combined text report.

    Parameters
    ----------
    final : List[str]
        length‐N list of bead‐tags (one per atom index)
    mapping : List[List[List[Any]]]
        full_mapping structure from your algorithm
    smiles : str
        the SMILES string used to position atoms
    compound_name : str
        base name for the output file (and written into the header)
    write_file : bool
        if True, writes out `{compound_name}.txt`; otherwise just returns data

    Returns
    -------
    beads_data : List[dict]
        list of per‐bead dicts with keys
         - index, bead_type, composed_of_index, composed_of_element,
           num_connections, mass, x,y,z
    connections_data : List[dict]
        list of per‐bond dicts with keys
         - index, connection_index, bond_type
    """
    # 1) group atoms into beads & get bead‐type labels
    coi, bead_types       = group_beads_by_type(final)
    # 2) get element lists per bead
    coe                   = get_composed_elements(coi, mapping)
    # 3) discover which beads are bonded to which
    origin, connected, bond_types = compute_beads_connections(coi, mapping)
    # 4) count total connections per bead
    num_conn              = compute_num_connections(len(coi), origin, connected)
    # 5) clean up bead‐type strings (strip +’s etc)
    fixed_types           = fix_beadtypes(bead_types)
    # 6) pull atomic coords from SMILES
    # --- TRY to get a 3D conformer; on failure, use (0,0,0) for all atoms ---
    try:
        atom_coords  = get_coordinates_from_smiles(smiles)
        coord_lookup = {
            c['index']: (c['x'] / 10.0, c['y'] / 10.0, c['z'] / 10.0)
            for c in atom_coords
        }
    except ValueError as e:
        print(f"[export_bead_mapping] warning: 3D embed failed ({e}), using zero‐coords")

        # Build a coord_lookup that maps every atom's global index to (0.0,0.0,0.0)
        coord_lookup = {}
        for sec in mapping:
            for atom_info in sec:
                gi = atom_info[0]   # global index
                coord_lookup[gi] = (0.0, 0.0, 0.0)
                
    # --- assemble beads_data ---
    beads_data = []
    for b_idx, atom_list in enumerate(coi):
        # centroid of its atoms
        xs, ys, zs = zip(*(coord_lookup[i] for i in atom_list))
        x0, y0, z0 = sum(xs)/len(xs), sum(ys)/len(ys), sum(zs)/len(zs)
        # mass by bead‐type prefix
        bt = fixed_types[b_idx]
        if bt.startswith('T'):
            mass = 36.0
        elif bt.startswith('S'):
            mass = 54.0
        else:
            mass = 72.0

        beads_data.append({
            'index':             b_idx,
            'bead_type':         bt,
            'composed_of_index': list(atom_list),
            'composed_of_element': list(coe[b_idx]),
            'num_connections':   num_conn[b_idx],
            'mass':              mass,
            'x':                  x0,
            'y':                  y0,
            'z':                  z0
        })

    # --- assemble connections_data ---
    connections_data = []
    for ai, aj, bond in zip(origin, connected, bond_types):
        connections_data.append({
            'index':            ai,
            'connection_index': aj,
            'bond_type':        bond
        })

    # --- write text report if requested ---
    if write_file:
        fmt_hdr = "{:<8} {:<10} {:<20} {:<20} {:<12} {:<8} {:<8} {:<8} {:<8}"
        fmt_ln  = "{:<8d} {:<10s} {:<20s} {:<20s} {:<12d} {:<8.1f} {:<8.4f} {:<8.4f} {:<8.4f}"
        lines = []
        lines.append(f"Compound: {compound_name}")
        lines.append(f"SMILES:   {smiles}\n")
        lines.append("All beads")
        lines.append("-"*80)
        lines.append(fmt_hdr.format(
            "index","bead_type","atoms_idx","atoms_elem",
            "#conns","mass","x","y","z"
        ))
        for b in beads_data:
            idxs = ",".join(str(i) for i in b['composed_of_index'])
            elems= ",".join(b['composed_of_element'])
            lines.append(fmt_ln.format(
                b['index'], b['bead_type'], idxs, elems,
                b['num_connections'], b['mass'],
                b['x'], b['y'], b['z']
            ))
        lines.append("\nBead connections")
        lines.append("-"*80)
        lines.append("{:<8} {:<8} {:<8}".format("index","to_idx","bond"))
        for c in connections_data:
            lines.append(f"{c['index']:<8d} {c['connection_index']:<8d} {c['bond_type']:<8}")

        with open(f"{compound_name}.txt","w") as f:
            f.write("\n".join(lines))

    return beads_data, connections_data
    
def gro_maker(beads_data, compound_name, box=None):
    """
    Generate a GROMACS .gro file for a single-residue system of beads.

    Parameters:
    - beads_data: list of dicts, each with keys:
        * index (int): zero-based bead index
        * bead_type (str): bead name/type
        * x, y, z (float): bead coordinates (in nm)
    - compound_name: str, used as residue name (max 5 chars recommended)
    - box: tuple of floats (bx, by, bz) for box vectors; if None, defaults to (5.0, 5.0, 5.0)

    Output:
    Writes a file named '{compound_name}.gro' with:
      Line 1: title = '1' + compound_name
      Line 2: number of beads
      Lines 3...: one line per bead:
        residue number (1), residue name, atom name, atom number, x, y, z
      Last line: box vectors
    """
    # ensure residue name fits in 5 chars
    res_name = "res"
    n_beads = len(beads_data)

    # default box if not given
    if box is None:
        box = (10.0, 10.0, 10.0)

    lines = []
    # Title
    lines.append(f"1{compound_name}")
    # Number of atoms
    lines.append(f"{n_beads}")
    
    # Atom lines
    count = 0
    for bead in beads_data:
        count = count + 1
        idx       = bead['index'] + 1           # GRO atom indices are 1-based
        #atom_name = bead['bead_type'][:5].rjust(5)
        x, y, z   = bead['x'], bead['y'], bead['z']
        # Format: %5d%5s%5s%5d%8.3f%8.3f%8.3f
        line = (
            f"{1:5d}{res_name:<5s}C{count}{idx:5d}"
            f"{x:8.3f}{y:8.3f}{z:8.3f}"
        )
        lines.append(line)

    # Box vector line: three floats %10.5f
    bx, by, bz = box
    lines.append(f"{bx:10.5f}{by:10.5f}{bz:10.5f}")

    # Write out
    filename = f"{compound_name}.gro"
    with open(filename, 'w') as f:
        f.write("\n".join(lines) + "\n")

    print(f"GRO file written to: {filename}")
    
def itp_maker(beads_data, connections_data, compound_name, write_file=True):
    """
    Generate a GROMACS .itp include file for a single-residue coarse-grained system.

    Parameters:
    - beads_data: list of dicts, each with keys:
        * index (int, 0-based)
        * bead_type (str)
        * mass (float)
    - connections_data: list of dicts, each with keys:
        * index (int, 0-based)
        * connection_index (int, 0-based)
        * bond_type (str or int)
    - compound_name: str, used as residue name and file basename
    - write_file: bool, if True writes '{compound_name}.itp' to disk

    Returns:
    - text: the full .itp content as a string
    """
    lines = []

    # ——————————————
    # Molecule type
    # ——————————————
    lines.append("[ moleculetype ]")
    lines.append("; name        nrexcl")
    lines.append("res           1")
    lines.append("")

    # ——————————————
    # Atoms
    # ——————————————
    lines.append("[ atoms ]")
    header_fmt = "; {:>5s} {:>5s} {:>5s} {:>5s} {:>5s} {:>5s} {:>10s} {:>10s}"
    lines.append(header_fmt.format("nr", "type", "resnr", "resid", "atom", "cgnr", "charge", "mass"))
    atom_fmt = "{:7d} {:5s} {:5d} {:5s} {:5s} {:5d} {:10.3f} {:10.3f}"
    count = 0
    for bead in beads_data:
        count = count + 1
        nr     = bead["index"] + 1
        atype  = bead["bead_type"][:5]
        resnr  = 1
        resid  = "res"
        atom   = f"C{count}"
        cgnr   = nr
        charge = 0.0
        mass   = bead["mass"]
        lines.append(atom_fmt.format(nr, atype, resnr, resid, atom, cgnr, charge, mass))
    lines.append("")

    # ——————————————
    # Bonds
    # ——————————————
    lines.append("[ bonds ]")
    lines.append(";  ai    aj   funct      b0(nm)        k")
    bond_fmt = "{:5d} {:5d} {:5d} {:12.5f} {:8d}"
    coord_lookup = { bead["index"]: (bead["x"], bead["y"], bead["z"])
                 for bead in beads_data }
    for conn in connections_data:
        ai    = conn["index"] + 1
        aj    = conn["connection_index"] + 1
        funct = 1
        # compute Euclidean distance b0 in nm
        x1, y1, z1 = coord_lookup[conn["index"]]
        x2, y2, z2 = coord_lookup[conn["connection_index"]]
        dx, dy, dz = x1 - x2, y1 - y2, z1 - z2
        b0 = math.sqrt(dx*dx + dy*dy + dz*dz)
          
        k = 20000  # constant force constant for all bonds
          
        lines.append(bond_fmt.format(ai, aj, funct, b0, k))
    lines.append("")

    text = "\n".join(lines)
    if write_file:
        fname = f"{compound_name}.itp"
        with open(fname, "w") as f:
            f.write(text)
        print(f"ITP file written to: {fname}")

    return text

