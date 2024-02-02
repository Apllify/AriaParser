import parser
import argparse

if __name__ == "__main__":
    argparser = argparse.ArgumentParser(description="A python tool for finding molecule conformations.")
    argparser.add_argument("prot_file", type = str, required  = True, metavar = ".prot file")
    argparser.add_argument("peaks_file", type = str, required  = True, metavar = ".peaks file")

    args = argparser.parse_args()

    #TODO : call parser routines to get usable values

    #TODO : use parser results to build .dat file

    #TODO : solve mp model with ampl and output .pdb file 