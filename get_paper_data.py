import numpy as np
from amplpy import AMPL, add_to_path
import aria_parser
import init_x
import utils.file_manip as file_manip
import os
import save_coordinates_pdb_format
import RMSD

def compare_init(par_file, top_file, prot_file, peaks_file, mod_file, solver):
    add_to_path(r'C:\Users\Admin\ampl_mswin64')
    #PAR PARSE
    with open(par_file, "r") as stream: 
        par_content = stream.read()

    cov_lengths, angle_lengths = aria_parser.parse_par(par_content)

    #TOP parse
    with open(top_file, "r") as stream:
        top_content = stream.read()
    res_info_dict = aria_parser.parse_top(top_content, cov_lengths, angle_lengths)

    #PEAKS PARSE
    with open(peaks_file, "r") as stream : 
        peaks_content = stream.read()

    #PROT PARSE
    with open(prot_file, "r") as stream : 
        prot_content = stream.read()
    _, _, res_id_to_AA_all = aria_parser.parse_prot(prot_content, res_info_dict)
    num_res_all = len(res_id_to_AA_all)

    for length in range(4, 11):
        average_iterations = [0, 0, 0]
        num_iter = 20
        for iter in range(num_iter):
            while 1:
                x = np.random.randint(num_res_all - length) + 1
                y = x + length - 1
                for i in range (x, y+1):
                    if res_id_to_AA_all[i] != 'XAA':
                        break
                else:
                    continue
                #if at least one residue is not XAA, we can proceed
                break
            terms = prot_file.split('.')
            new_prot_file = f'{terms[0]}_{x}_to_{y}.{terms[1]}'  # Output file path
            file_manip.extract_residues_in_range(prot_file, new_prot_file, x, y)
            #PROT PARSE
            with open(new_prot_file, "r") as stream : 
                prot_content = stream.read()
            os.remove(new_prot_file)
            chem_shift_to_atom, atom_set, res_id_to_AA = aria_parser.parse_prot(prot_content, res_info_dict)

            #store NOE measurements with the associated atoms
            noe_assign = aria_parser.parse_peaks(peaks_content, chem_shift_to_atom)

            #now get all non-NOE distances
            cov_assign, ang_assign = aria_parser.compute_dists(atom_set, res_info_dict, res_id_to_AA)

            for init in range(3): # 0: no init, 1: init normal, 2: init multiple
                #initialize x values
                if init == 0 or init == 1:
                    atoms_to_initial_coord = init_x.initialize_x(atom_set, res_id_to_AA)
                elif init == 2:
                    atoms_to_initial_coord = init_x.initialize_x_multiple_aas(atom_set, res_id_to_AA)
                    
                #calculate the scale factor
                # sf = calc_sf.calc_sf(atom_set, noe_assign, res_id_to_AA)
                #write out our results
                aria_parser.write_data(atom_set, noe_assign, cov_assign, ang_assign, atoms_to_initial_coord, "model/NOE_generated.dat")

                ampl = AMPL()
                ampl.read(mod_file)
                ampl.read_data('model/NOE_generated.dat')
                if init == 0:
                    ampl.eval("for {a in ATOMS, d in DIM} { let x[a, d] := Uniform(-M, M); }")
                elif init >= 1:
                    ampl.eval("let {a in ATOMS, d in DIM} x[a, d] := init_x[a, d]; # explicit initialization")
                ampl.eval(f'option solver {solver};')
                # ampl.eval(f'option {solver}_options "max_iter=6000";')  # Replace 100 with your desired maximum number of iterations')
                solve_output = ampl.get_output("solve;")
                objectives_length = 0
                for line in solve_output.split('\n'):
                    terms = line.split()
                    if len(terms) != 0 and (terms[0].isnumeric() or (terms[0][-1] == 'r' and terms[0][:-1].isnumeric())) :
                        # objectives.append(terms[1])
                        objectives_length += 1
                    elif len(terms) != 0 and terms[0] == 'EXIT:':
                        print(f'length = {length}, iter = {iter}, init = {init}, ans = {terms[1]}')
                # print(f'For init type: {init}, length: {length}, number of iteration = {len(objectives)}')
                # print(init, objectives_length)
                average_iterations[init] += objectives_length

                with open(f"output_{init}_{length}", "w+") as f:
                    f.write(solve_output)

        average_iterations = list(map(lambda x: x / num_iter, average_iterations))
        print(f'For length {length}, average iterations: {average_iterations}')


def compare_RMSD(par_file, top_file, prot_file, peaks_file, first_mod, second_mod, solver):
    add_to_path(r'C:\Users\Admin\ampl_mswin64')
    #PAR PARSE
    with open(par_file, "r") as stream: 
        par_content = stream.read()

    cov_lengths, angle_lengths = aria_parser.parse_par(par_content)

    #TOP parse
    with open(top_file, "r") as stream:
        top_content = stream.read()
    res_info_dict = aria_parser.parse_top(top_content, cov_lengths, angle_lengths)

    #PEAKS PARSE
    with open(peaks_file, "r") as stream : 
        peaks_content = stream.read()

    #PROT PARSE
    with open(prot_file, "r") as stream : 
        prot_content = stream.read()
    _, _, res_id_to_AA_all = aria_parser.parse_prot(prot_content, res_info_dict)
    num_res_all = len(res_id_to_AA_all)

    for length in range(1, 31):
        average_RMSD = [0, 0]
        num_iter = 10
        for iter in range(num_iter):
            while 1:
                x = np.random.randint(num_res_all - length) + 1
                y = x + length - 1
                for i in range (x, y+1):
                    if res_id_to_AA_all[i] != 'XAA':
                        break
                else:
                    continue
                #if at least one residue is not XAA, we can proceed
                break
            terms = prot_file.split('.')
            new_prot_file = f'{terms[0]}_{x}_to_{y}.{terms[1]}'  # Output file path
            file_manip.extract_residues_in_range(prot_file, new_prot_file, x, y)
            #PROT PARSE
            with open(new_prot_file, "r") as stream : 
                prot_content = stream.read()
            chem_shift_to_atom, atom_set, res_id_to_AA = aria_parser.parse_prot(prot_content, res_info_dict)

            #store NOE measurements with the associated atoms
            noe_assign = aria_parser.parse_peaks(peaks_content, chem_shift_to_atom)

            #now get all non-NOE distances
            cov_assign, ang_assign = aria_parser.compute_dists(atom_set, res_info_dict, res_id_to_AA)
            atoms_to_initial_coord = init_x.initialize_x(atom_set, res_id_to_AA)
            
            aria_parser.write_data(atom_set, noe_assign, cov_assign, ang_assign, atoms_to_initial_coord, "model/NOE_generated.dat")

            ampl = AMPL()
            ampl.read(first_mod)
            ampl.read_data('model/NOE_generated.dat')
            # ampl.eval("for {a in ATOMS, d in DIM} { let x[a, d] := Uniform(-M, M); }")
            ampl.eval("let {a in ATOMS, d in DIM} x[a, d] := init_x[a, d]; # explicit initialization")
            ampl.eval(f'option solver {solver};')
            first_solve_output = ampl.get_output("solve;")
            ampl.eval(f'display x > first_model_output.txt;')
            save_coordinates_pdb_format.parse_output('first_model_output.txt', par_file, top_file, new_prot_file)
            # normalized model
            ampl = AMPL()
            ampl.read(second_mod)
            ampl.read_data('model/NOE_generated.dat')
            # ampl.eval("for {a in ATOMS, d in DIM} { let x[a, d] := Uniform(-M, M); }")
            ampl.eval("let {a in ATOMS, d in DIM} x[a, d] := init_x[a, d]; # explicit initialization")
            ampl.eval(f'option solver {solver};')
            second_solve_output = ampl.get_output("solve;")
            ampl.eval(f'display x > second_model_output.txt;')
            save_coordinates_pdb_format.parse_output('second_model_output.txt', par_file, top_file, new_prot_file)
            average_RMSD[0] += RMSD.RMSDcalc_pdb('first_model_output.pdb')
            average_RMSD[1] += RMSD.RMSDcalc_pdb('second_model_output.pdb')
            os.remove(new_prot_file)
        average_RMSD = list(map(lambda x: x / num_iter, average_RMSD))
        print(f'length:{length}, average_RMSD: {average_RMSD}')

if '__main__':
    # compare_init('data/aria.par', 'data/aria.top', 'data/hmqcnoe.prot', 'data/hmqcnoe.peaks', 'model/NOE_normalize.mod', 'ipopt')
    # file_manip.extract_residues_in_range('data/hmqcnoe.prot', 'data/hmqcnoe_10_19.prot', 10, 19)
    compare_RMSD('data/aria.par', 'data/aria.top', 'data/hmqcnoe.prot', 'data/hmqcnoe.peaks', 'model/NOE.mod', 'model/NOE_normalize.mod', 'ipopt')