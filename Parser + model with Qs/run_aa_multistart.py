import aria_parser
import init_x
import evaluate_model as ev
from amplpy import AMPL, Environment, DataFrame, add_to_path
import math
import pandas as pd
import save_coordinates_pdb_format as sv
import numpy as np

"""Program that runs amino acid random restart and collects data. 
To be used together with init_x and save_coordinates_pdb_format from AriaParser
Rk. To use, change names of files where appropriate """

numitr = 25
estar = math.inf
xstar = None
sdiststar = None
smstar = None
spstar = None
fail = 0


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

errors = []
filenames = ["data/hmqcnoe_33_37.prot", "data/hmqcnoe_45_49.prot", "data/hmqcnoe_60_5.prot", "data/hmqcnoe_85_5.prot",
             "data/hsqcnoe_33_37.prot", "data/hsqcnoe_45_49.prot", "data/hsqcnoe_60_5.prot", "data/hsqcnoe_85_5.prot"]
n = len(filenames)
peaks = ["m", "m", "m", "m", "s", "s", "s", "s"]

for i in range(n):
    print(f'processing {i} aas')
    estar = math.inf
    Xstar = None
    fail = 0
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

    atoms_to_initial_coord = init_x.initialize_x(atom_set, res_id_to_AA)

    #write out our results
    aria_parser.write_data(atom_set, noe_assign, cov_assign, ang_assign, q_dists, "model/aa_ms_inpt.dat")

    

    for _ in range(numitr):
        if fail == 5:
            break

        atoms_to_initial_coord = init_x.initialize_x(atom_set, res_id_to_AA)
        aria_parser.write_initx(atoms_to_initial_coord)

    
        ampl.read("/Users/anjamatic/AriaParser/model/NOE_aams.run")
        E = ampl.get_value("estar")

        
        if E == 0:
            fail += 1
            continue
    

        if E < estar:
            estar = E
            xstar = ampl.get_data("xstar").toPandas()
            sdiststar = ampl.get_data("sdist").to_pandas()
            smstar = ampl.get_data("sm").to_pandas()
            spstar = ampl.get_data("sp").to_pandas()

    errors.append(estar)
    xstar.to_csv("/Users/anjamatic/AriaParser/xstar_otpt.csv")

    with open(f'xstar_otpt.csv', "r") as stream: 
                model_output = stream.read()
    atom_names, atom_ress, Xs = sv.parse_output_csv(model_output)
    AAs = [res_id_to_AA[res] for res in atom_ress]

    #sort according to residue number
    zipped = zip(AAs, atom_ress, atom_names, Xs)
    zipped = sorted(zipped, key=lambda pair: pair[1])
    AAs, atom_ress, atom_names, Xs = map(list, zip(*zipped))
    sv.save_coordinates_pdb_format(f'PDB_structs/{i}_aarrs_{file[5:]}.pdb', "PDB_ID", "A", AAs, atom_ress, atom_names, np.array(Xs), "METHOD")

    sdiststar.to_csv(f'/Users/anjamatic/AriaParser/output_data/sdists_aarrs_{file[5:]}')
    smstar.to_csv(f'/Users/anjamatic/AriaParser/output_data/sm_aarrs_{file[5:]}')
    spstar.to_csv(f'/Users/anjamatic/AriaParser/output_data/sp_aarrs_{file[5:]}')
    
print(errors)
df = pd.DataFrame(errors)
df.to_csv("/Users/anjamatic/AriaParser/errors_aars_5_AA.csv")


