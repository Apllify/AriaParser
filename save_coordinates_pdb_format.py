#!/usr/bin/python3

from datetime import date
import numpy as np

## ------------------------------------------------------------------------------------------------------
# From Wagner's code
def save_coordinates_pdb_format(filename, pdb_id, chain_name, res_name, res_seq, atom_name, X, method):

   n = len(X)

   with open(filename, 'w') as new_file:
      new_file.write('HEADER                                                {}   {}         \n'.format(date.today(), pdb_id))

      if method != 'PDB':
         new_file.write('TITLE     {} solution                                                          \n'.format(method))
      else:
         new_file.write('TITLE     PDB structure                                                         \n')

      new_file.write('REMARK   2                                                                      \n')
      new_file.write('REMARK   2 RESOLUTION.    2.00 ANGSTROMS.                                       \n')
      new_file.write('MODEL        1                                                                  \n')

      for k in range(n):
         atomType   = "ATOM  "
         serial     = k+1
         name       = atom_name[k]
         altLoc     = " "
         resName    = res_name[k]
         chainID    = chain_name
         resSeq     = int(res_seq[k])
         iCode      = " "
         x          = X[k,0]
         y          = X[k,1]
         z          = X[k,2]
         occupancy  = 1
         tempFactor = 1
         element    = atom_name[k][0]
         charge     = "  "

         new_file.write('{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s}   {:8.3f}{:8.3f}{:8.3f}{:6.2f}{:6.2f}          {:>2s}{:2s}\n'.format(atomType, serial, name, altLoc, resName, chainID, resSeq, iCode, x, y, z, occupancy, tempFactor, element, charge))
         #new_file.write('ATOM   {:4.0f} {} {} {}{:4.0f}    {:8.3f}{:8.3f}{:8.3f} {:5.2f} {:5.2f}           {}  \n'.format('ATOM', k+1, a_name, res_name[k], chain_name, res_seq[k], X[k,0], X[k,1], X[k,2], 1, 1, atom_name[k][0]))

      atomType = "TER   "
      serial   = k+1
      name     = "    "
      altLoc   = " "
      resName  = res_name[k]
      chainID  = chain_name
      resSeq   = int(res_seq[k])
      iCode    = " "

      new_file.write('{:6s}{:5d} {:4s}{:1s}{:3s} {:1s}{:4d}{:1s}\n'.format(atomType, serial, name, altLoc, resName, chainID, resSeq, iCode))
      #new_file.write('TER    {:4.0f} {} {} {}{:4.0f}\n'.format(k+1, '    ', res_name[k], chain_name, res_seq[k]))

def parse_output(content : str) -> tuple[list, list, list]:
   """
   Transforms an AMPL model output into a PDB file
   """
   lines = content.split("\n")
   atom_names = []
   atom_ress = []
   Xs = []

   for line in lines : 
      terms = line.split()
      #skip invalid lines
      if line == "" or line[0] in ['!', ':', ';'] or len(terms) != 4:
         continue

      #get all info from current line
      try : 
         atom = terms[0]
         x, y, z = float(terms[1]), float(terms[2]), float(terms[3])
         Xs.append([x, y, z])
      except : 
         print("failure reading data")
         continue

      #remember this atom's name and residue for later
      atom_terms = atom.split('_')
      atom_name, atom_res = atom_terms[0].upper(), atom_terms[1]
      atom_names.append(atom_name)
      atom_ress.append(atom_res)

   return (atom_names, atom_ress, Xs)

with open("model_output.txt", "r") as stream: 
   model_output = stream.read()
atom_names, atom_ress, Xs = parse_output(model_output)
# Place holder name
save_coordinates_pdb_format("output.pdb", "PDB_ID", "CHAIN_NAME", atom_ress, [-1 for _ in range(len(atom_names))], atom_names, np.array(Xs), "METHOD")