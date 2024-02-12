import aria_parser
import argparse

from pprint import pprint # debugging

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="A python tool for finding molecule conformations.")
    argparser.add_argument("prot_file", type = str, metavar = ".prot file")
    argparser.add_argument("peaks_file", type = str, metavar = ".peaks file")
    argparser.add_argument("par_file", type = str, metavar = ".par file")
    argparser.add_argument("top_file", type = str, metavar = ".top file")

    args = argparser.parse_args()

    #PAR PARSE
    with open(args.par_file, "r") as stream : 
        par_content = stream.read()
   

    d1_assign, d2_assign = aria_parser.parse_par(par_content)
    #pprint(d1_assign)
    #pprint(d2_assign)

     
    #TOP parse
    with open(args.top_file, "r") as stream:
        top_content = stream.read()
    res_info_dict = aria_parser.parse_top(top_content, d1_assign)


    #PROT PARSE
    with open(args.prot_file, "r") as stream : 
        prot_content = stream.read()
    chem_shift_to_atom, atom_set = aria_parser.parse_prot(prot_content, res_info_dict)
    #print(f'ASP: {res_info_dict["ASP"]}')


    #PEAKS PARSE
    with open(args.peaks_file, "r") as stream : 
        peaks_content = stream.read()

    #store NOE measurements with the associated atoms
    noe_assign = aria_parser.parse_peaks(peaks_content, chem_shift_to_atom)
    #print(noe_assign)

    cov_dists = aria_parser.compute_dists(atom_set, dict())
    ##pprint(cov_dists)

    
    #write out our results
    #aria_parser.write_data(atom_set, noe_assign, cov_dists, "test_file.txt")




    #TODO : use parser results to build .dat file

    #TODO : solve mp model with ampl and output .pdb file 