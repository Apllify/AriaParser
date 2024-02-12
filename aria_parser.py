import math
from collections import OrderedDict

#Type defs
Atom = str #  ID of the form ATOMTYPE_RESNUM

ChemShiftToAtom = dict[tuple[float, float], Atom]
AtomsByRes = dict[int, set[Atom]] #atoms by residue

PairAssignment = dict[ tuple[Atom, Atom], float ]
TripletAssignment = list[ tuple[ list[Atom], list[Atom], list[Atom], float ] ]


ResInfo = tuple[set[Atom],   PairAssignment] 
#              (res_content, res_bond_lengths)
ResInfoDict = dict[str, ResInfo]




#GLOBAL DO NOT TOUCH
CACHE_SIZE = 30
lookup_cache = OrderedDict()


def parse_prot(content : str, res_info_dict : ResInfoDict) -> tuple[ChemShiftToAtom, AtomsByRes]:
    """
    Transforms a protein info file in to a dictionary
    which maps each chem shift interval : (LB, UB) to 
    a protein id 
    """

    chem_shift_to_atom = dict()
    atom_set = dict()

    lines = content.split("\n")
    last_res_id = 1

    seq = set() #atom sequence of current residue
    AA_name = ""
    hydrogen_counter = 0 #unify notations

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

            atom_type = terms[3]
            res_id = int(terms[4])
        except : 
            continue

        #remember this atom's name and ID for later
        atom_name =  f"{atom_type.upper()}_{res_id}"
        if atom_set.get(res_id) : 
            atom_set[res_id].add(atom_name)
        else : 
            atom_set[res_id] = {atom_name}

        #also store the chem shift of atom if it is well defined
        if chem_shift != 999 : 
            chem_shift_to_atom[(chem_shift - err, chem_shift + err)] = atom_name


        #check if we just switched to a new residue
        if last_res_id != res_id:
            AA_name = match_AA(seq, res_info_dict)
            if AA_name != "???": 
                print(AA_name)
            else: 
                print(f'res_nr: {last_res_id}, AA not found')
                print(seq)
            
            atom_set[res_id].add(f'O_{res_id}')
            last_res_id = res_id
            seq = set()
            hydrogen_counter = 0
        
        if "Q" in atom_type:
            hydrogen_counter = 0
            group_id = ""
            continue
       
        elif atom_type[-1].isdigit() and atom_type[0] == "H" and f'C{atom_type[1:]}' not in seq and f'N{atom_type[1:]}' not in seq:
            hydrogen_counter += 1
            seq.add(f'{atom_type[:-1]}{hydrogen_counter}')
            continue
        else:
            hydrogen_counter = 0
            seq.add(atom_type)
            continue
        
        
    return (chem_shift_to_atom, atom_set)

def atoms_from_shift(chem_shift : float, chem_shift_to_atom : ChemShiftToAtom) -> list[Atom] : 
    """
    Determine the id of an atom from its chemical shift
    using the value intervals established from the prot file 

    Returns the empty string when no valid atom is found
    """
    global lookup_cache

    atoms = []

    #check if we already computed this recently
    if chem_shift in lookup_cache:
        lookup_cache.move_to_end(chem_shift)
        return lookup_cache[chem_shift]

    #iterate through all intervals
    for (lb, ub) in chem_shift_to_atom.keys():
        if chem_shift >= lb and chem_shift <= ub: 
            atoms.append(chem_shift_to_atom[lb, ub])
    

    #store result and return
    if len(lookup_cache) >= CACHE_SIZE:
        lookup_cache.popitem(last=False)

    lookup_cache[chem_shift] = atoms
    
    return atoms

def match_AA(atoms: set[Atom], res_info_dict: ResInfoDict) -> str:
    """
    Tries to find the amino acid which most closely
    matches the given atom set.
    
    Returns "???" if no match found
    """
    for name, (seq, _) in res_info_dict.items():
        if seq == atoms:
           #special case since CYS and SER have same composition
            if name in ("CYS", "SER"):
               return "CYS/SER" 
           
            return name
    
    if atoms | {"HZ3"} == res_info_dict["LYS"][0]: return "LYS"
    if atoms | {"HE2"} == res_info_dict["HIS"][0]: return "HIS"
    if atoms - {"HE1"} == res_info_dict["GLU"][0]: return "GLU"
    if atoms - {"HD1"} == res_info_dict["ASP"][0]: return "ASP"
    if {"NE", "NH1", "NH2"} <= atoms: return "ARG"

    return "???"
    
def parse_peaks(content : str, chem_shift_to_atom : ChemShiftToAtom) -> TripletAssignment:
    """
    Uses the prot information to read the peaks file 
    and transform each triplet of chem shifts into their 
    respective atom IDs and stores the NOE measure for each.
    """
    
    lines = content.split("\n")
    #Assignmnet of rho to corresponding atoms follows scheme detailed in README.txt
    triplet_assignment : TripletAssignment = []
    triplets = []
    unbound_carbon_triplet = []
    hydrogen_pairs = dict()

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
            m1, m2, m3 = map(int, terms[10:13]) #atoms
        except : 
            continue

        #ignore negative peaks
        if noe <= 0 : 
            continue

        #ignore lines that don't have 3 unique atoms / have a zero atom
        uniques = set([chem1, chem2, chem3])
        if len(uniques) <= 2 or 0 in uniques : 
            continue 
        
        #find the atom IDs from their chem shifts
        Ip = atoms_from_shift(chem1, chem_shift_to_atom)
        Iq = atoms_from_shift(chem2, chem_shift_to_atom)
        Ir = atoms_from_shift(chem3, chem_shift_to_atom)

        #make sure only C or Ns are in Ip
        Ip = list(filter(lambda atom : atom[0] in "CN", Ip))

        #make sure only Hs and Qs are in Ip,Ir
        Iq = list(filter(lambda atom : atom[0] in "HQ", Iq))
        Ir = list(filter(lambda atom : atom[0] in "HQ", Ir))
        
        #Check that each chemical shift correspnds to something
        if 0 in (len(Ip), len(Iq), len(Ir)):
            continue
    

        triplet_assignment.append((Ip, Iq, Ir, noe))

    return triplet_assignment

def parse_par(content : str) -> tuple[PairAssignment, PairAssignment]: 
    """
    Parses the content of the .par file
    into a list of generic distances between atom types

    Returns a tuple where elem 1 is the assignment obtained from 
    bond statements, while elem2 is the assignment obtained
    from angle statements.
    """

    bond_assignment : PairAssignment = dict()
    angle_assignment : PairAssignment = dict()

    for line in content.split("\n") : 

        #file is case insensitive
        line = line.upper()
        terms = line.split()
        # remove terms of {sd=...} in aria.par
        terms = list(filter(lambda term: ('{' not in term) and ('}' not in term), terms))

        if len(terms) == 0 : 
            continue

        opcode = terms[0]
        match opcode :  
            case "BOND":
                #store the length of the given bond
                try : 
                    m1, m2 = terms[1], terms[2]
                    dist = float(terms[4])
                except : 
                    continue
                bond_assignment[(m1, m2)] = dist

            case "ANGLE" : 
                #retrieve angle data
                try : 
                    m1, m2, m3 = terms[1], terms[2], terms[3]
                    angle = float(terms[5])
                except : 
                    continue

                #get the two lengths of the triangle (we check both orders)
                d1 = bond_assignment.get((m1, m2)) or bond_assignment.get((m2, m1))
                d2 = bond_assignment.get((m1, m3)) or bond_assignment.get((m3, m1))
                if d1 == None or d2 == None : 
                    continue

                #compute and store the 3rd length
                dist = math.sqrt( d1**2 + d2**2 - d1*d2*math.cos(math.radians(angle)) ) 
                angle_assignment[(m2, m3)] = dist
    return (bond_assignment, angle_assignment)



def parse_top(content: str, bond_lengths: PairAssignment) -> ResInfoDict :
    """
    Parses the top file and finds the sequence and
    bond information for each amino acid.

    The word residue is used to refer to AAs for
    consistency and simplicity. 
    """
    
    res_info_dict : ResInfoDict = dict()

    inside_res = False
    res_name = ""
    res_atoms = set()
    res_bond_lengths : PairAssignment = dict()

    res_atom_to_type = dict() #doesn't appear in output


    for line in content.split('\n'):
        line = line.upper()
        terms = line.split()
        if len(terms) == 0:
            continue

        opcode = terms[0]

        if opcode == "RESIDUE" : 
            inside_res = True
            res_name = terms[1]

            if res_name == "ACE": 
                break #rest of file is irrelevant 



        #only process commands within a residue
        if inside_res : 
            match opcode:
               
                case "ATOM":
                    name, type= terms[1], terms[2]

                    if name[0] not in ("O", "S"): 
                        res_atoms.add(name)
                    res_atom_to_type[name] = type[5:]
                
                case "BOND":
                    n = len(terms)
                    for i in range(0, n, 3):
                        #ignore comments
                        if terms[i][0] == "!": 
                            continue

                        #find the length of bond from .par data
                        a1, a2 = terms[i+1], terms[i+2]
                        general_bond = (res_atom_to_type[a1], res_atom_to_type[a2])
                        general_bond = tuple(sorted(general_bond))

                        #not all cov lengths are known
                        length = bond_lengths.get(general_bond)
                        if length != None :
                            res_bond_lengths[(a1, a2)] = length

                case "END":
                    #add our new entry and reset variables
                    res_info_dict[res_name] = (res_atoms, res_bond_lengths)

                    inside_res = False
                    res_name = ""
                    res_atoms = set()
                    res_bond_lengths = dict()

                    res_atom_to_type = dict()


    return res_info_dict
           

def compute_dists(atoms : AtomsByRes, generic_dists : PairAssignment) -> PairAssignment :
    """
    Uses the generic covalent bond distances 
    obtained from the .par file to find all the
    distances applicable to our atoms set
    """
    
    #for now only use some hardcoded generic dists
    backbone_dists = {
        ("n", "hn") : 0.980,
        ("n", "ca") : 1.458,
        ("ca", "ha") : 1.080,
        ("ca", "c") : 1.525,
        ("c", "o") : 1.231 #oxygen on carbon backbone, does not appear in .prot file
    }

    atom_dists : PairAssignment = dict()

    #systematically does not att C_id and O
    for res_id, residue in atoms.items():
        for (a1, a2) in backbone_dists.keys():
            #check for each known distance pair if it is in this residue
            a1_spec = f"{a1}_{res_id}"
            a2_spec = f"{a2}_{res_id}"

            if a1_spec in residue and a2_spec in residue: 
                atom_dists[(a1_spec, a2_spec)] = backbone_dists[(a1, a2)]
            
    return atom_dists

def write_data(atoms: AtomsByRes, rhos: TripletAssignment,  cov_dists: PairAssignment, filename = "NOE_data.dat"):

    with open(filename, "w") as outfile:
        #define atoms set
        outfile.write("set ATOMS := ")
        residues = atoms.values()
        for residue in residues : 
            for a in residue: outfile.write(f' {a}')
        outfile.write(";\n")

        #define covalent distances set
        outfile.write("set COVDISTS := ")
        for a1, a2 in cov_dists.keys(): outfile.write(f' {a1} {a2}')
        outfile.write(";\n")

        #define distance set
        outfile.write("set DISTS := ")
        for _, a1, a2, _ in rhos: outfile.write(f' {a1} {a2}')
        outfile.write(";\n")

        #define RHOS set
        outfile.write("set RHOS := ")
        for a1, a2, a3, _ in rhos: outfile.write(f' {a1} {a2} {a3}')
        outfile.write(";\n")

        #give covalent distance data
        outfile.write("param  CovDists := ")
        for (a1, a2), dist in cov_dists.items(): outfile.write(f' {a1} {a2} {dist}')
        outfile.write(";\n")

        #give RHO data
        outfile.write("param rho := ")
        for a1, a2, a3, chemshift in rhos: outfile.write(f' {a1} {a2} {a3} {chemshift}')
        outfile.write(";\n")




            




