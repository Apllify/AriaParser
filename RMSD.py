import sys
import Bio.PDB as PDB
from os import listdir
from os.path import isfile, join
import numpy as np
import argparse

# Store atom name and 3d coordinates
class Atom:
    def __init__(self, name = "", coord = np.array(3)):
        self.name = name
        self.coord = coord
    def __lt__(self, other):
        return self.name < other.name
    def __repr__(self):
        return f"Atom({self.name}, {self.coord})"

def get_aa_names_from_file(file):
    """
    Get all aa names in a top file
    """
    with open(f"{file}", "r") as stream:
        content = stream.read()
    aa_names = []
    res_name = ""

    for line in content.split('\n'):
        line = line.upper()
        terms = line.split()
        if len(terms) == 0:
            continue
        opcode = terms[0]
        if opcode == "RESIDUE": 
            res_name = terms[1]
            if res_name == "ACE": 
                break #rest of file is irrelevant 
            if res_name not in ['ZN', 'CHEX', 'HISH']:
                aa_names.append(res_name)
    return aa_names
    
def RMSDcalc(X0, Xr):
    """
    Calculate the RMSD between two numpy arrays
    """
    n = X0.shape[0]
    X0bar = np.zeros((n, 3))
    Xrbar = np.zeros((n, 3))
    cmX0 = np.sum(X0, axis=0)/n
    cmXr = np.sum(Xr, axis=0)/n

    for i in range(n):
        X0bar[i, :] = X0[i, :] - cmX0
        Xrbar[i, :] = Xr[i, :] - cmXr

    U, S, Vt = np.linalg.svd(np.dot(np.transpose(Xrbar), X0bar))
    Q = np.dot(np.transpose(Vt), np.transpose(U))
    Xrbar = np.dot(Xrbar, np.transpose(Q))
    RMSD = np.linalg.norm(X0bar - Xrbar, ord='fro') / np.sqrt(n)

    return X0bar, Xrbar, RMSD

def RMSDcalc_res(res0, resr):
    """
    Calculate the RMSD between res0 (true residue) and resr (generated residue)
    """
    res0_atoms = [Atom(atom.name, atom.coord) for atom in res0]
    resr_atoms = [Atom(atom.name, atom.coord) for atom in resr]

    res0_atoms = sorted(res0_atoms)
    for i in range(len(res0_atoms)):
        # Change all 2-3 hydrogen atoms to 1-2 hydrogen atoms
        if i < len(res0_atoms)-2 \
        and res0_atoms[i+1].name[0] == 'H' and res0_atoms[i+1].name[-1] == '2' \
        and res0_atoms[i+2].name[0] == 'H' and res0_atoms[i+2].name[-1] == '3' \
        and res0_atoms[i].name != res0_atoms[i+1].name[:-1] + '1':
            res0_atoms[i+1].name = res0_atoms[i+1].name[:-1] + '1'
            res0_atoms[i+2].name = res0_atoms[i+2].name[:-1] + '2'
        # Change all H atoms to HN atoms
        if res0_atoms[i].name == 'H':
            res0_atoms[i].name = 'HN'
    res0_atoms = sorted(res0_atoms)

    X0 = np.array([x.coord for x in res0_atoms])
    Xr = np.array([x.coord for x in resr_atoms])
    return RMSDcalc(X0, Xr)

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description="Calculating the RMSD between a generated pdb and residue position from pdb folder.")
    argparser.add_argument("gen_pdb", type = str, metavar = "Generated pdb")
    argparser.add_argument("num_res", type = int, metavar= "Number of residues to compare with")

    args = argparser.parse_args()

    #TOP parse
    aa_names = get_aa_names_from_file("data/aria.top")
    residues = {x: [] for x in aa_names}
    count_aa_found = 0

    # Parse downloaded pdb
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
                    if len(residues[res_name]) < args.num_res:
                        count_aa_found += 1
                        residues[res_name].append(residue)
        if count_aa_found == len(aa_names) * args.num_res:
            break

    # Parse generated pdb
    struct = pdbparser.get_structure(args.gen_pdb.split('.')[0], args.gen_pdb)
    for model in struct:
        for chain in model:
            for res_id, residue in enumerate(chain):
                gen_res_name = residue.get_resname()
                RMSD_avg = 0
                if gen_res_name not in residues:
                    print(f'For amino acid {gen_res_name} at residue {res_id+1}, not found in any pdb file')
                    continue

                for true_residue in residues[gen_res_name]:
                    X0bar, Xrbar, RMSD = RMSDcalc_res(true_residue, residue)
                    RMSD_avg += RMSD
                RMSD_avg /= len(residues[gen_res_name])
                print(f'For amino acid {gen_res_name} at residue {res_id+1}, average RMSD = {RMSD_avg}')
