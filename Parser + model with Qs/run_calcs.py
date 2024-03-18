import aria_parser
import evaluate_model as ev
from amplpy import AMPL, Environment, DataFrame, add_to_path
import math
import numpy as np
import save_coordinates_pdb_format as sv
import pandas as pd

errors = []
Xstar = None
fail = 0
"""Program to run NOE models and collect data.
Rk: change paths were appropriate"""

#PAR PARSE
with open("/Users/anjamatic/AriaParser/data/aria.par", "r") as stream : 
    par_content = stream.read()
   
cov_lengths, angle_lengths = aria_parser.parse_par(par_content)

#TOP parse
with open("/Users/anjamatic/AriaParser/data/aria.top", "r") as stream:
    top_content = stream.read()
res_info_dict = aria_parser.parse_top(top_content, cov_lengths, angle_lengths)



add_to_path(r"/Users/anjamatic/ampl") 
ampl = AMPL()

filenames = ["data/hmqcnoe_33_37.prot", "data/hmqcnoe_45_49.prot", "data/hmqcnoe_60_5.prot", "data/hmqcnoe_85_5.prot",
             "data/hsqcnoe_33_37.prot", "data/hsqcnoe_45_49.prot", "data/hsqcnoe_60_5.prot", "data/hsqcnoe_85_5.prot"]
n = len(filenames)
peaks = ["m", "m", "m", "m", "s", "s", "s", "s"]

for mod in ["", "_vns","_ms"]:
    errors = []
    fail = 0
    for i in range(n):
        estar = math.inf
        file = filenames[i]
    
        #PROT PARSE
        with open(file, "r") as stream : 
            prot_content = stream.read()
        chem_shift_to_atom, atom_set, res_id_to_AA, res_id_to_Q = aria_parser.parse_prot(prot_content, res_info_dict)

        #PEAKS PARSE
        with open(f'/Users/anjamatic/AriaParser/data/h{peaks[i]}qcnoe.peaks', "r") as stream : 
            peaks_content = stream.read()

        #store NOE measurements with the associated atoms
        noe_assign = aria_parser.parse_peaks(peaks_content, chem_shift_to_atom)

        #now get all non-NOE distances
        cov_assign, ang_assign = aria_parser.compute_dists(atom_set, res_info_dict, res_id_to_AA)

        #and q distances
        q_dists = aria_parser.compute_Q(res_id_to_Q, atom_set, res_id_to_AA, res_info_dict)

        #write out our results
        aria_parser.write_data(atom_set, noe_assign, cov_assign, ang_assign, q_dists, "model/aa_ms_inpt.dat")
        
        ampl.read(f'/Users/anjamatic/AriaParser/model/NOE{mod}.run')
        E = ampl.get_value("estar")
        errors.append(E)

        Xstar = ampl.get_data("x").to_pandas()
        Xstar.to_csv("/Users/anjamatic/AriaParser/xstar_otpt.csv")

        ampl.get_data("sdist").to_pandas().to_csv(f'/Users/anjamatic/AriaParser/output_data/sdists{mod}_{file[5:]}')
        ampl.get_data("sm").to_pandas().to_csv(f'/Users/anjamatic/AriaParser/output_data/sm{mod}_{file[5:]}')
        ampl.get_data("sp").to_pandas().to_csv(f'/Users/anjamatic/AriaParser/output_data/sp{mod}_{file[5:]}')

        with open(f'xstar_otpt.csv', "r") as stream: 
                model_output = stream.read()
        atom_names, atom_ress, Xs = sv.parse_output_csv(model_output)
        AAs = [res_id_to_AA[res] for res in atom_ress]

        #sort according to residue number
        zipped = zip(AAs, atom_ress, atom_names, Xs)
        zipped = sorted(zipped, key=lambda pair: pair[1])
        AAs, atom_ress, atom_names, Xs = map(list, zip(*zipped))
        sv.save_coordinates_pdb_format(f'PDB_structs/{i}{mod}_{file[5:]}.pdb', "PDB_ID", "A", AAs, atom_ress, atom_names, np.array(Xs), "METHOD")

    df = pd.DataFrame(errors)
    df.to_csv(f'/Users/anjamatic/AriaParser/output_data/errors{mod}_5_AA.csv')


