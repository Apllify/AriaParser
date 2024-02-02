import parser
import argparse

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="A python tool for finding molecule conformations.")
    argparser.add_argument("prot_file", type = str, metavar = ".prot file")
    argparser.add_argument("peaks_file", type = str, metavar = ".peaks file")

    args = argparser.parse_args()

    with open(args.prot_file, "r") as stream : 
        prot_content = stream.read()
    chem_shift_to_id, id_to_mol = parser.parse_prot(prot_content)


    with open(args.peaks_file, "r") as stream : 
        peaks_content = stream.read()

    p_assign = parser.parse_peaks(peaks_content, chem_shift_to_id, id_to_mol)
    print(str(p_assign).replace(",", "\n"))

    #TODO : use parser results to build .dat file

    #TODO : solve mp model with ampl and output .pdb file 