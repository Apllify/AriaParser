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
    # print(f'N new position: {rotated_positions[N_idx]}, C new position: {rotated_positions[C_idx]}, V: {V}')
    
    return rotated_positions

# Example usage:
# V = np.array([x, y, z])  # Your target vector
# rotated_positions = rotate_molecule_to_vector(atom_positions, N_idx, C_idx, V)
# print(rotated_positions)

def initialize_x(atom_set, res_id_to_AA, dim=3):
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
    
    # Find atom to coord mapping
    residues_atom_set = list(atom_set.values())
    maxCovDist = 2 # from aria.par
    bound = sum([len(residue) for residue in residues_atom_set]) * maxCovDist / 2 / 10
    atoms_to_coord = {atom: np.random.uniform(low=-bound, high=bound, size=dim) for residue in residues_atom_set for atom in residue}
    visited = dict()
    cur_pos = np.array([0,0,0])
    for res_id in atom_set:
        aa = res_id_to_AA[res_id]
        if aa not in residues:
            # already initialized to uniformly random values so we can skip
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

        while 1:
            # let each residue occupies a random cube in space
            # such that consecutive residue is adjacent
            random_V =  np.random.choice([-1,0,1], size=3, replace=True)
            if tuple(cur_pos + random_V) not in visited:
                visited[tuple(cur_pos + random_V)] = 1
                cur_pos += random_V
                break

        atom_coords = np.array([atom.coord for atom in residue_atoms])
        rotated_atom_coords = rotate_molecule_to_vector(atom_coords, N_idx, C_idx, random_V)
        center = np.mean(rotated_atom_coords, axis=0)
        # first normalization, so that each residue's center is at the origin point
        normalize_atom_coords = rotated_atom_coords - center

        RES_SIZE = 5 # changeable
        for i, atom in enumerate(residue_atoms):
            atoms_to_coord[f'{atom.name}_{res_id}'] = normalize_atom_coords[i] + cur_pos * RES_SIZE

    # final normalization
    mean = np.mean(np.array(list(atoms_to_coord.values())), axis = 0)
    atoms_to_coord = {key: val - mean for key, val in atoms_to_coord.items()}
    return atoms_to_coord