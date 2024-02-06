import math
from collections import OrderedDict

#Type defs
Atom = str #  ID of the form ATOMTYPE_RESNUM

ChemShiftToAtom = dict[tuple[float, float], Atom]
AtomSet = set[Atom]

PairAssignment = dict[ tuple[Atom, Atom], float ]
TripletAssignment = list[ tuple[ Atom, Atom, Atom, float ] ]


#GLOBAL DO NOT TOUCH
CACHE_SIZE = 30
lookup_cache = OrderedDict()


def parse_prot(content : str) -> tuple[ChemShiftToAtom, AtomSet]:
    """
    Transforms a protein info file in to a dictionary
    which maps each chem shift interval : (LB, UB) to 
    a protein id 
    """

    chem_shift_to_mol = dict()
    mol_set = set()

    lines = content.split("\n")

    for line in lines : 
        terms = line.split()

        #skip invalid lines
        if line == "" or line[0] == "!" or len(terms) != 5 : 
            continue


        #get all info from current line
        try : 
            id = int(terms[0])
            chem_shift = float(terms[1])
            err =  float(terms[2])

            mol_type = terms[3]
            res_id = int(terms[4])
        except : 
            continue

        #remember this molecule's name and ID for later
        mol_name =  f"{mol_type}_{res_id}"
        mol_set.add(mol_name)

        #also store the chem shift of molecule if it is well defined
        if chem_shift != 999 : 
            chem_shift_to_mol[(chem_shift - err, chem_shift + err)] = mol_name

        
    return (chem_shift_to_mol, mol_set)

        


def mol_from_shift(chem_shift : float, chem_shift_to_mol : ChemShiftToAtom) -> Atom : 
    """
    Determine the id of an atom from its chemical shift
    using the value intervals established from the prot file 

    Returns the empty string when no valid mol is found
    """
    global lookup_cache

    cur_min_name = ""
    cur_min =  math.inf

    #check if we already computed this recently
    if chem_shift in lookup_cache:
        lookup_cache.move_to_end(chem_shift)
        return lookup_cache[chem_shift]

    #iterate through all intervals
    for (lb, ub) in chem_shift_to_mol.keys():
        if chem_shift >= lb and chem_shift <= ub: 

            dist = abs(chem_shift - (ub+lb)/2)
            if  dist < cur_min:
                cur_min_name = chem_shift_to_mol[(lb, ub)]
                cur_min = dist

    #store result and return
    if len(lookup_cache) >= CACHE_SIZE:
        lookup_cache.popitem(last=False)

    lookup_cache[chem_shift] = cur_min_name
    
    return cur_min_name



def parse_peaks(content : str, chem_shift_to_mol : ChemShiftToAtom) -> TripletAssignment:
    """
    Uses the prot information to read the peaks file 
    and transform each triplet of chem shifts into their 
    respective atom IDs and stores the NOE measure for each.
    """
    
    lines = content.split("\n")
    #Assignmnet of rho to corresponding atoms follows scheme detailed in README.txt
    triplet_assignment : TripletAssignment = []

    for line in lines : 
        #ignore comments
        if line.strip() == "" or line.strip()[0] == "#":
            continue

        #get info from line
        try : 
            terms = line.split()
            id = int(terms[0])
            chem1, chem2, chem3 = map(float, terms[1:4]) #chem shifts
            noe = float(terms[6])
            m1, m2, m3 = map(int, terms[10:13]) #molecules
        except : 
            continue

        #ignore negative peaks
        if noe <= 0 : 
            continue

        #ignore lines that don't have 3 unique molecules / have a zero mol
        uniques = set([chem1, chem2, chem3])
        if len(uniques) <= 2 or 0 in uniques : 
            continue 
        
        #find the atom IDs from their chem shifts
        Ip = mol_from_shift(chem1, chem_shift_to_mol)
        Iq = mol_from_shift(chem2, chem_shift_to_mol)
        Ir = mol_from_shift(chem3, chem_shift_to_mol)
        
        if -1 in (Ip, Iq, Ir):
            continue

        #store the result
        triplet_assignment.append( (Ip, Iq, Ir, noe) )

    return triplet_assignment



def parse_par(content : str) -> PairAssignment: 
    """
    Parses the content of the .par file
    into a list of distances
    """

    pair_assignment : PairAssignment = dict()

    for line in content : 

        #file is case insensitive
        line = line.lower()
        terms = line.split()

        opcode = terms[0]

        match opcode :  

            case "bond":
                #store the length of the given bond
                try : 
                    m1, m2 = terms[1], terms[2]
                    dist = float(terms[6])
                except : 
                    continue

                pair_assignment[(m1, m2)] = dist

            case "angle" : 
                #retrieve angle data
                try : 
                    m1, m2, m3 = terms[1], terms[2], terms[3]
                    angle = float(terms[7])
                except : 
                    continue

                #get the two lengths of the triangle (we check both orders)
                d1 = pair_assignment.get((m1, m2)) or pair_assignment.get((m2, m1))
                d2 = pair_assignment.get((m1, m3)) or pair_assignment.get((m3, m1))
                if d1 == None or d2 == None : 
                    continue

                #compute and store the 3rd length
                dist = math.sqrt( d1**2 + d2**2 - d1*d2*math.cos(math.radians(angle)) ) 
                pair_assignment[(m2, m3)] = dist

    return pair_assignment

def write_data(atoms: AtomSet, rhos: TripletAssignment, filename = "NOE_data.dat"):

    with open(filename, "w") as outfile:
        #define atoms set
        outfile.write("set ATOMS := ")
        outfile.write(*atoms)
        outfile.write(";\n")

        #define covalent distances set
        outfile.write("set COVDISTS := ;")

        #define distance set
        outfile.write("set DISTS := ")
        for a1, a2, _, _ in rhos: outfile.write(f'{a1} {a2}')
        outfile.write(";\n")

        #define RHOS set
        outfile.write("set RHOS := ")
        for a1, a2, a3, _ in rhos: outfile.write(f'{a1} {a2} {a3}')
        outfile.write(";\n")

        #give RHO data
        outfile.write("param rho := ")
        for a1, a2, a3, chemshift in rhos: outfile.write(f'{a1} {a2} {a3} {chemshift}')
        outfile.write(";\n")




            




