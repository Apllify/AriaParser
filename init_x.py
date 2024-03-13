import Bio.PDB as PDB
from os import listdir
from os.path import isfile, join
import numpy as np
import utils.atom_class as atom_class
import utils.atom_func as atom_func

RES_SIZE = 5 # changeable
# ASSUME ALL CYS/SER ARE CYS
three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'XAA' : 'C'}
one_to_three = {value: key for (key, value) in three_to_one.items()}

def rotate_molecule_to_vector(atom_positions, N_idx, C_idx, V):
    # Translate the molecule so that N is at the origin
    N = atom_positions[N_idx]
    translated_positions = atom_positions - N
    
    # Compute the NC vector and normalize it
    NC_vector = translated_positions[C_idx]
    NC_norm = NC_vector / np.linalg.norm(NC_vector)
    
    # Normalize the target vector V
    V_norm = V / np.linalg.norm(V)
    
    # Compute the rotation axis (cross product between NC and V)
    rotation_axis = np.cross(NC_norm, V_norm)
    rotation_axis_norm = rotation_axis / np.linalg.norm(rotation_axis)
    
    # Compute the angle between NC and V
    cos_angle = np.dot(NC_norm, V_norm)
    angle = np.arccos(cos_angle)
    
    # Rodrigues' rotation formula components
    K = np.array([[0, -rotation_axis_norm[2], rotation_axis_norm[1]],
                  [rotation_axis_norm[2], 0, -rotation_axis_norm[0]],
                  [-rotation_axis_norm[1], rotation_axis_norm[0], 0]])
    I = np.eye(3)
    
    # Compute the rotation matrix
    rotation_matrix = I + np.sin(angle) * K + (1 - np.cos(angle)) * (K @ K)
    
    # Apply the rotation matrix to all atom positions
    rotated_positions = translated_positions @ rotation_matrix.T
    
    return rotated_positions

def initialize_x(atom_set, res_id_to_AA, dim=3):
    """
    Initialize positions according to amino acid internal positions from pdb files
    """
    aa_names = list(set(res_id_to_AA.values()))
    if 'XAA' in aa_names:
        aa_names.remove('XAA')
    residues = {x: [] for x in aa_names}
    count_aa_found = 0

    # Parse downloaded pdb
    num_res = 1
    path = "data/pdb"
    files = [f for f in listdir(path) if isfile(join(path, f))]
    for file in files:
        pdbparser = PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
        struct = pdbparser.get_structure(file.split('.')[0], f"data/pdb/{file}")
        for model in struct:
            for chain in model:
                for residue in chain:
                    res_name = residue.get_resname()
                    if res_name not in residues:
                        continue
                    if len(residues[res_name]) < num_res:
                        count_aa_found += 1
                        residues[res_name].append(residue)
        if count_aa_found == len(aa_names) * num_res:
            break
    
    # Find atom to coordinate mapping
    residues_atom_set = list(atom_set.values())
    maxCovDist = 2 # from aria.par
    bound = sum([len(residue) for residue in residues_atom_set]) * maxCovDist / 2 / 60
    atoms_to_coord = {atom: np.random.uniform(low=-bound, high=bound, size=dim) for residue in residues_atom_set for atom in residue}
    cur_pos = np.array([0, 0, 0])
    visited = {(0, 0, 0): 1}
    for res_id in atom_set:
        aa = res_id_to_AA[res_id]
        # Already initialized to uniformly random values so we can skip
        if aa not in residues:
            continue
        residue_atoms = atom_func.get_correct_atom_name(residues[aa][0], res_id, atoms_to_coord)

        # Find C and N to align the molecule's C-N vector along a random vector
        C_idx = -1
        N_idx = -1
        for i, atom in enumerate(residue_atoms):
            if atom.name == f'C_{res_id}':
                C_idx = i
            elif atom.name == f'N_{res_id}':
                N_idx = i

        # Let each residue occupies a random cube in space
        # such that consecutive residue is adjacent
        while 1:
            random_V =  np.random.choice([-1,0,1], size=3, replace=True)
            if tuple(cur_pos + random_V) not in visited:
                visited[tuple(cur_pos + random_V)] = 1
                cur_pos += random_V
                break

        atom_coords = np.array([atom.coord for atom in residue_atoms])
        rotated_atom_coords = rotate_molecule_to_vector(atom_coords, N_idx, C_idx, random_V)
        center = np.mean(rotated_atom_coords, axis=0)
        # First normalization, so that each residue's center is at the origin point
        normalized_atom_coords = rotated_atom_coords - center

        for i, atom in enumerate(residue_atoms):
            atoms_to_coord[atom.name] = normalized_atom_coords[i] + cur_pos * RES_SIZE

    # Final normalization
    mean = np.mean(np.array(list(atoms_to_coord.values())), axis = 0)
    atoms_to_coord = {key: val - mean for key, val in atoms_to_coord.items()}
    return atoms_to_coord

def get_condensed_pdb_file(res_id_to_AA):
    all_pdb_seqs = []
    path = "data/pdb"
    files = [f for f in listdir(path) if isfile(join(path, f))]
    for file in files:
        pdbparser = PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
        struct = pdbparser.get_structure(file.split('.')[0], f"data/pdb/{file}")
        pdb_seq = [file, ' ']
        for model in struct:
            if model.id != 1:
                continue
            for chain in model:
                # Only consider A chain
                if chain.id != 'A':
                    continue
                for residue in chain.get_residues():
                    if residue.get_resname() in three_to_one:
                        pdb_seq.append(three_to_one[residue.get_resname()])
        pdb_seq = ''.join(pdb_seq)
        all_pdb_seqs.append(pdb_seq)
    all_pdb_seqs = '\n'.join(all_pdb_seqs)
    file = open('data/pdb_condensed.txt','w+')
    file.write(all_pdb_seqs)

def initialize_x_multiple_aas(atom_set, res_id_to_AA, dim = 3):
    # Only generate file if not existing
    if not isfile("data/pdb_condensed.txt"):
        get_condensed_pdb_file(res_id_to_AA)

    prot_ress = list(res_id_to_AA.values())
    prot_ress_one = [three_to_one[x] for x in prot_ress]
    s_prot =  ''.join(prot_ress_one)

    with open('data/pdb_condensed.txt') as f:
        lines = f.readlines()
    files_seqs = []
    for line in lines:
        terms = line.split()
        if len(terms) < 2:
            continue
        files_seqs.append((terms[0], terms[1]))
    pdbparser = PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
    
    # Find atom to coordinate mapping
    residues_atom_set = list(atom_set.values())
    maxCovDist = 2 # from aria.par
    bound = sum([len(residue) for residue in residues_atom_set]) * maxCovDist / 2 / 60
    atoms_to_coord = {atom: np.random.uniform(low=-bound, high=bound, size=dim) for residue in residues_atom_set for atom in residue}
    cur_pos = np.array([0, 0, 0])
    visited = {(0, 0, 0): 1}

    i = 0
    while i < len(s_prot):
        file_with_seq = ""
        pos = 0
        for j in range(i+1, len(s_prot)+2):
            for file_seq in files_seqs:
                file, seq = file_seq
                if seq.find(s_prot[i:j]) != -1:
                    pos = seq.find(s_prot[i:j])
                    file_with_seq = file
                    break
            else:
                break
        struct = pdbparser.get_structure(file.split('.')[0], f"data/pdb/{file_with_seq}")
        for model in struct:
            if model.id != 1:
                continue
            for chain in model:
                if chain.id != 'A':
                    continue
                local_res_group = []
                for res_id, residue in enumerate(chain):
                    prot_res_id = res_id - pos + i + 1
                    if prot_res_id >= i + 1 and prot_res_id < j:
                        residue_atoms = atom_func.get_correct_atom_name(residue, prot_res_id, atoms_to_coord)
                        local_res_group += residue_atoms

                atom_coords = np.array([atom.coord for atom in local_res_group])

                # Find C and N to align the molecule's C-N vector along a random vector
                C_idx = -1
                N_idx = -1
                for id, atom in enumerate(local_res_group):
                    if atom.name == f'C_{i+1}':
                        C_idx = id
                    elif atom.name == f'N_{j-1}':
                        N_idx = id

                # Let each local residue group occupies a random cube in space
                # such that consecutive residue is adjacent
                while 1:
                    random_V =  np.random.choice([-1,0,1], size=3, replace=True)
                    if tuple(cur_pos + random_V) not in visited:
                        visited[tuple(cur_pos + random_V)] = 1
                        cur_pos += random_V
                        break

                atom_coords = np.array([atom.coord for atom in local_res_group])
                rotated_atom_coords = rotate_molecule_to_vector(atom_coords, N_idx, C_idx, random_V)
                center = np.mean(rotated_atom_coords, axis=0)
                # First normalization, so that each residue's center is at the origin point
                normalized_atom_coords = rotated_atom_coords - center
  
                for i, atom in enumerate(local_res_group):
                    atoms_to_coord[atom.name] = normalized_atom_coords[i] + cur_pos * RES_SIZE
        i = j-1

    # Final normalization
    mean = np.mean(np.array(list(atoms_to_coord.values())), axis = 0)
    atoms_to_coord = {key: val - mean for key, val in atoms_to_coord.items()}
    return atoms_to_coord

    # # brute force all Cys/Ser
    # for i in range(1<<6):
    #     for k in range(6):
    #         assert(s_prot[cys_ser_id[k] - 1] in ['X', 'C', 'S'])
    #         if i & (1<<k):
    #             s_prot = s_prot[:cys_ser_id[k]-1] + 'C' + s_prot[cys_ser_id[k]:]
    #         else:
    #             s_prot = s_prot[:cys_ser_id[k]-1] + 'S' + s_prot[cys_ser_id[k]:]
    #         # print(s_prot)

    # max_length = max(max_length, pylcs.lcs_string_length(s_prot, s_seq))
    # print(pylcs.lcs_string_length(s_prot, s_seq), s_seq)
    # break