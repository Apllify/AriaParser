import sys
import Bio.PDB as PDB
import Bio
from os import listdir
from os.path import isfile, join
import numpy as np
import RMSD

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
    #TOP parse
    aa_names = RMSD.get_aa_names_from_file("data/aria.top")
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
    visited = dict()
    cur_pos = np.array([0,0,0])
    for res_id in atom_set:
        aa = res_id_to_AA[res_id]
        # Already initialized to uniformly random values so we can skip
        if aa not in residues:
            continue
        residue = residues[aa][0]
        residue_atoms = [RMSD.Atom(atom.name, atom.coord) for atom in residue]
        residue_atoms = sorted(residue_atoms)
        for i in range(len(residue_atoms)):
            # Change all 2-3 hydrogen atoms to 1-2 hydrogen atoms
            if i < len(residue_atoms)-2 \
            and residue_atoms[i+1].name[0] == 'H' and residue_atoms[i+1].name[-1] == '2' \
            and residue_atoms[i+2].name[0] == 'H' and residue_atoms[i+2].name[-1] == '3' \
            and residue_atoms[i].name != residue_atoms[i+1].name[:-1] + '1':
                residue_atoms[i+1].name = residue_atoms[i+1].name[:-1] + '1'
                residue_atoms[i+2].name = residue_atoms[i+2].name[:-1] + '2'
            # Change all H atoms to HN atoms
            if residue_atoms[i].name == 'H':
                residue_atoms[i].name = 'HN'
        residue_atoms = sorted([atom for atom in residue_atoms if f'{atom.name}_{res_id}' in atoms_to_coord])

        # Find C and N to align the molecule's C-N vector along a random vector
        C_idx = -1
        N_idx = -1
        for i, atom in enumerate(residue_atoms):
            if atom.name == 'C':
                C_idx = i
            elif atom.name == 'N':
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

        RES_SIZE = 5 # changeable
        for i, atom in enumerate(residue_atoms):
            atoms_to_coord[f'{atom.name}_{res_id}'] = normalized_atom_coords[i] + cur_pos * RES_SIZE

    # Final normalization
    mean = np.mean(np.array(list(atoms_to_coord.values())), axis = 0)
    atoms_to_coord = {key: val - mean for key, val in atoms_to_coord.items()}
    return atoms_to_coord

def initialize_x_multiple_aas(atom_set, res_id_to_AA, dim = 3):
    path = "data/pdb"
    files = [f for f in listdir(path) if isfile(join(path, f))]
    prot_ress = list(res_id_to_AA.values())
    for i, x in enumerate(prot_ress):
        if x == 'XAA':
            print(i+1)
    three_to_one = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
     'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
     'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
     'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M', 'XAA' : 'C'}
    one_to_three = {value: key for (key, value) in three_to_one.items()}
    max_length = 0
    prot_ress_one = [three_to_one[x] for x in prot_ress]
    s_prot =  ''.join(prot_ress_one)
    print(s_prot)
    cys_ser_id = [x+1 for x in range(len(s_prot)) if s_prot[x] == 'X']
    all_pdb_seqs = []
    for file in files:
        pdbparser = PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
        struct = pdbparser.get_structure(file.split('.')[0], f"data/pdb/{file}")
        for model in struct:
            for chain in model:
                print(chain.get_residues())
        # read_seqres = 0
        # seqress = []
        # with open(f'{path}/{file}') as f:
        #     while 1:
        #         line = f.readline()
        #         if not read_seqres and line[:6] == 'SEQRES':
        #             read_seqres = 1
        #         elif read_seqres and line[:6] != 'SEQRES':
        #             read_seqres = 0
        #             break
        #         if read_seqres:
        #             terms = line.split()
        #             # ignore sequence that is not A
        #             if terms[2] != 'A':
        #                 continue
        #             seqress += [x for x in terms[4:] if x[0] != 'D' and len(x) == 3]
        # seqress_one = [three_to_one[x] if x in three_to_one else '0' for x in seqress]
        # s_seq = ''.join (seqress_one)
        # all_pdb_seqs.append(s_seq)
    
    # all_pdb_seqs = '\t'.join(all_pdb_seqs)
    # print(len(all_pdb_seqs))
    # print(len(s_prot))
    # i = 0
    # while 1:
    #     broken = 0
    #     seq_file = ""
    #     for j in range(i+1, len(s_prot)+2):
    #         for file_num, seq in enumerate(all_pdb_seqs):
    #             if s_prot[i:j] in seq:
    #                 seq_file = files[file_num]
    #                 break
    #         else:
    #             break
    #     print(i, j-1, seq_file)
    #     i = j-1
    #     if i == len(s_prot):
    #         break

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
    print(max_length)