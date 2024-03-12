import Bio.PDB as PDB
from os import listdir
from os.path import isfile, join
import utils.atom_func as atom_func
import numpy as np

def calc_sf(atom_set, noe_assign, res_id_to_AA, num_res = 1):
    all_res = set(res_id_to_AA.values())
    residues = {x: [] for x in all_res}
    count_aa_found = 0
    path = "data/pdb"
    files = [f for f in listdir(path) if isfile(join(path, f))]
    pdbparser = PDB.PDBParser(QUIET=True)   # suppress PDBConstructionWarning
    for file in files:
        struct = pdbparser.get_structure(file.split('.')[0], f"data/pdb/{file}")
        for model in struct:
            for chain in model:
                for residue in chain:
                    res_name = residue.get_resname()
                    if res_name not in all_res:
                        continue
                    if len(residues[res_name]) < num_res:
                        count_aa_found += 1
                        residues[res_name].append(atom_func.get_correct_atom_name_no_filter(residue, 0))
        if count_aa_found == len(all_res) * num_res:
            break
    res_inner_dist = dict()
    for res_name, residues_one_type in residues.items():
        res_inner_dist[res_name] = dict()
        for residue in residues_one_type:
            for atom1 in residue:
                if atom1.name[0] != 'H':
                    continue
                for atom2 in residue:
                    if atom2.name[0] != 'H':
                        continue
                    if atom1 == atom2:
                        continue
                    res_inner_dist[res_name][atom1.name.split('_')[0] + atom2.name.split('_')[0]] = np.linalg.norm(atom1.coord - atom2.coord)

    avg_sf = 0
    cnt = 0
    for (Ip, Iq, Ir, noe) in noe_assign:
        noe_calc = 0

        for q in Iq:
            res_id_q = q.split('_')[1]
            for r in Ir:
                if q == r:
                    continue
                combined_atom_name = q.split('_')[0] + r.split('_')[0]
                res_id_r = r.split('_')[1]
                if res_id_q != res_id_r:
                    continue
                aa = res_id_to_AA[int(res_id_q)]
                if combined_atom_name not in res_inner_dist[aa]:
                    continue
                noe_calc += 1/(res_inner_dist[aa][combined_atom_name] ** 6)
        if noe_calc == 0:
            continue
        avg_sf += noe/noe_calc
        cnt += 1
    avg_sf /= cnt
    print(avg_sf)