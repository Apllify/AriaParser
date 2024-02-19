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
    
# Calculate the RMSD
def RMSDcalc(X0, Xr):
   n = X0.shape[0]
   X0bar = np.zeros((n,3))
   Xrbar = np.zeros((n,3))
   cmX0 = np.sum(X0, axis=0)/n
   cmXr = np.sum(Xr, axis=0)/n

   for i in range(n):
      X0bar[i,:] = X0[i,:] - cmX0
      Xrbar[i,:] = Xr[i,:] - cmXr

   U,S,Vt = np.linalg.svd(np.dot(np.transpose(Xrbar), X0bar))
   Q = np.dot(np.transpose(Vt), np.transpose(U))
   Xrbar = np.dot(Xrbar, np.transpose(Q))
   RMSD = np.linalg.norm(X0bar - Xrbar, ord='fro')/np.sqrt(n)

   return X0bar, Xrbar, RMSD

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description="Calculating the RMSD between a generated pdb and residue position from pdb folder.")
    argparser.add_argument("gen_pdb", type = str, metavar = "generated pdb")

    args = argparser.parse_args()

    #TOP parse
    aa_names = get_aa_names_from_file("data/aria.top")
    res_atoms = {x: [] for x in aa_names}
    aa_found = 0

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
                    if res_name not in res_atoms:
                        continue
                    if len(res_atoms[res_name]) == 0:
                        aa_found += 1
                        res_atoms[res_name] = [Atom(atom.name, atom.coord) for atom in residue]
        if aa_found == len(aa_names):
            break

    # Parse generated pdb
    struct = pdbparser.get_structure(file.split('.')[0], args.gen_pdb)
    for model in struct:
        for chain in model:
            for residue in chain:
                gen_res_name = residue.get_resname()
                gen_atoms = sorted(Atom(atom.name, atom.coord) for atom in residue)

                true_atoms = sorted(res_atoms[gen_res_name])
                for i in range(len(true_atoms)):
                    # Change all 2-3 hydrogen atoms to 1-2 hydrogen atoms
                    if i < len(true_atoms)-2 \
                    and true_atoms[i+1].name[0] == 'H' and true_atoms[i+1].name[-1] == '2' \
                    and true_atoms[i+2].name[0] == 'H' and true_atoms[i+2].name[-1] == '3' \
                    and true_atoms[i].name != true_atoms[i+1].name[:-1] + '1':
                        true_atoms[i+1].name = true_atoms[i+1].name[:-1] + '1'
                        true_atoms[i+2].name = true_atoms[i+2].name[:-1] + '2'
                    if true_atoms[i].name == 'H':
                        true_atoms[i].name = 'HN'
                true_atoms = sorted(res_atoms[gen_res_name])
                X0 = np.array([x.coord for x in true_atoms])
                Xr = np.array([x.coord for x in gen_atoms])
                X0bar, Xrbar, RMSD = RMSDcalc(X0, Xr)
                print(f'RMSD: {RMSD}')
