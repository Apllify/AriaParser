import Bio.PDB as PDB
from os import listdir
from os.path import isfile, join
import numpy as np
import argparse
import utils.atom_class as atom_class
import utils.atom_func as atom_func
    
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
    res0_atoms = [atom_class.Atom(atom.name, atom.coord) for atom in res0]
    resr_atoms = [atom_class.Atom(atom.name, atom.coord) for atom in resr]

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

    if res0.get_resname() == 'ARG':
        res0_atoms = [x for x in res0_atoms if x.name != 'HH12']
    if res0.get_resname() == 'LYS':
        res0_atoms = [x for x in res0_atoms if x.name != 'HZ3']
    if res0.get_resname() == 'HIS':
        res0_atoms = [x for x in res0_atoms if x.name != 'HE2']
    if resr.get_resname() == 'GLU':
        resr_atoms = [x for x in resr_atoms if x.name != 'HE1']
    if resr.get_resname() == 'ASP':
        resr_atoms = [x for x in resr_atoms if x.name != 'HD1']
    # print(len(res0_atoms), len(resr_atoms))
    if len(res0_atoms) != len(resr_atoms):
        return "diff length", "diff length", "diff length"

    X0 = np.array([x.coord for x in res0_atoms])
    Xr = np.array([x.coord for x in resr_atoms])
    return RMSDcalc(X0, Xr)


def RMSDcalc_pdb(gen_pdb, num_res=10):
    """
    Calculating the RMSD between a generated pdb and random residue (that matches aa) from pdb folder
    """
    #TOP parse
    aa_names = atom_func.get_aa_names_from_file("data/aria.top")
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
                    if len(residues[res_name]) < num_res:
                        count_aa_found += 1
                        residues[res_name].append(residue)
        if count_aa_found == len(aa_names) * num_res:
            break

    # Parse generated pdb
    struct = pdbparser.get_structure(gen_pdb.split('.')[0], gen_pdb)
    sum_RMSD_avg = 0
    for model in struct:
        for chain in model:
            for res_id, residue in enumerate(chain):
                gen_res_name = residue.get_resname()
                RMSD_avg = 0
                if gen_res_name not in residues:
                    # print(f'For amino acid {gen_res_name} at residue {res_id+1}, not found in any pdb file')
                    continue

                cnt = 0
                for true_residue in residues[gen_res_name]:
                    # first/last residue have one more atom
                    # print(len(true_residue), len(residue), true_residue, residue)
                    X0bar, Xrbar, RMSD = RMSDcalc_res(true_residue, residue)
                    if RMSD != "diff length":
                        RMSD_avg += RMSD
                        cnt += 1
                RMSD_avg /= cnt
                sum_RMSD_avg += RMSD_avg
                # print(f'For amino acid {gen_res_name} at residue {res_id+1}, average RMSD = {RMSD_avg}')
    return sum_RMSD_avg

if __name__ == '__main__':
    argparser = argparse.ArgumentParser(description="Calculating the RMSD between a generated pdb and residue position from pdb folder.")
    argparser.add_argument("gen_pdb", type = str, metavar = "Generated pdb")
    argparser.add_argument("num_res", type = int, metavar= "Number of residues to compare with")

    args = argparser.parse_args()
    print(RMSDcalc_pdb(args.gen_pdb, args.num_res))
