import math
from collections import OrderedDict


#Type defs
ChemShiftToID = dict[tuple[float, float], int]
IDToMol = dict[int, str]

Mol = str | int # either ID or name
TripletAssignment = list[ tuple[ Mol, Mol, Mol, float ] ]


#GLOBAL DO NOT TOUCH
CACHE_SIZE = 30
lookup_cache = OrderedDict()


def parse_prot(content : str) -> tuple[ChemShiftToID, IDToMol]:
    """
    Transforms a protein info file in to a dictionary
    which maps each chem shift interval : (LB, UB) to 
    a protein id 
    """

    chem_shift_to_id = dict()
    id_to_mol = dict()

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
        id_to_mol[id] = mol_name

        #also store the chem shift of molecule if it is well defined
        if chem_shift != 999 : 
            chem_shift_to_id[(chem_shift - err, chem_shift + err)] = id

        
    return (chem_shift_to_id, id_to_mol)

        


def id_from_shift(chem_shift : float, chem_shift_to_id : ChemShiftToID) -> int : 
    """
    Determine the id of an atom from its chemical shift
    using the value intervals established from the prot file 

    Returns -1 when no valid id is found
    """
    global lookup_cache

    cur_min_id = -1
    cur_min =  math.inf

    #check if we already computed this recently
    if chem_shift in lookup_cache:
        lookup_cache.move_to_end(chem_shift)
        return lookup_cache[chem_shift]

    #iterate through all intervals
    for (lb, ub) in chem_shift_to_id.keys():
        if chem_shift >= lb and chem_shift <= ub: 

            dist = abs(chem_shift - (ub+lb)/2)
            if  dist < cur_min:
                cur_min_id = chem_shift_to_id[(lb, ub)]
                cur_min = dist

    #store result and return
    if len(lookup_cache) >= CACHE_SIZE:
        lookup_cache.popitem(last=False)

    lookup_cache[chem_shift] = cur_min_id
    
    return cur_min_id



def parse_peaks(content : str, chem_shift_to_id : ChemShiftToID, id_to_mol : IDToMol) -> TripletAssignment:
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
        Ip = id_from_shift(chem1, chem_shift_to_id)
        Iq = id_from_shift(chem2, chem_shift_to_id)
        Ir = id_from_shift(chem3, chem_shift_to_id)
        
        if -1 in (Ip, Iq, Ir):
            continue

        #convert the ids to names
        Np = id_to_mol.get(Ip)
        Nq = id_to_mol.get(Iq)
        Nr = id_to_mol.get(Ir)

        #store the result
        triplet_assignment.append( (Np, Nq, Nr, noe) )

    return triplet_assignment
