import sys
import Bio.PDB as PDB
from os import listdir
from os.path import isfile, join
import numpy as np
import RMSD

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
        residue_atoms = sorted(residue_atoms)
        for atom in residue_atoms:
            if f'{atom.name}_{res_id}' in atoms_to_coord:
                atoms_to_coord[f'{atom.name}_{res_id}'] = atom.coord

    return atoms_to_coord