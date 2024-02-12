import math
from collections import OrderedDict

#Type defs
Atom = str #  ID of the form ATOMTYPE_RESNUM

ChemShiftToAtom = dict[tuple[float, float], Atom]
AtomSet = dict[int, set[Atom]] #atoms by residue

PairAssignment = dict[ tuple[Atom, Atom], float ]
TripletAssignment = list[ tuple[ list, list, list, float ] ]


#GLOBAL DO NOT TOUCH
CACHE_SIZE = 30
lookup_cache = OrderedDict()


def parse_prot(content : str, AA_dict, AA_to_seq) -> tuple[ChemShiftToAtom, AtomSet]:
    """
    Transforms a protein info file in to a dictionary
    which maps each chem shift interval : (LB, UB) to 
    a protein id 
    """

    chem_shift_to_atom = dict()
    atom_set = dict()

    lines = content.split("\n")
    last_res_id = 1

    seq = set() #atom sequence of current resideu
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


        if last_res_id != res_id:
            
            AA_name = match_AA(seq, AA_to_seq)
            if AA_name != -1: print(AA_name)
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

def atom_from_shift(chem_shift : float, chem_shift_to_atom : ChemShiftToAtom) -> set() : 
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

def match_AA(target_AA: str, AA_to_seq: dict):
   
    for name, seq in AA_to_seq.items():
       if seq == target_AA:
           if name == "CYS" or name == "SER":
               return "CYS/SER"
           return name
    
    if target_AA |{"HZ3"} == AA_to_seq["LYS"]: return "LYS"
    if target_AA | {"HE2"} == AA_to_seq["HIS"]: return "HIS"
    if target_AA - {"HE1"} == AA_to_seq["GLU"]: return "GLU"
    if target_AA - {"HD1"} == AA_to_seq["ASP"]: return "ASP"
    if {"NE", "NH1", "NH2"} <= target_AA: return "ARG"

    return -1
    
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
        Ip = atom_from_shift(chem1, chem_shift_to_atom)
        Iq = atom_from_shift(chem2, chem_shift_to_atom)
        Ir = atom_from_shift(chem3, chem_shift_to_atom)
        
        #Check that each chemical shift correspnds to something
        if len(Ip) == 0 or len(Iq) == 0 or len(Ir) == 0:
            continue
        
        #To Rali and Cong: this was the best way I could think to do this lol 
        #Check that Ip only contains carbons or nitrogens
        if len([a for a in Ip if a[0] not in "CN"]) != 0:
            continue

        #Check that only hydrogens are in h and q
        if len([a for a in [*Iq, *Ir] if a[0] not in "HQ"]) != 0:
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

def parse_top(content: str, bond_lenghts: PairAssignment):
    """Parses the top file and finds the sequence + bond information for each amino acid"""
    seq_to_AA = dict()
    AA_name = ""
    AA_to_seq = dict()
    seq_set = set()
    seq = ""
    name_to_type = dict()
    bonds = []

    for line in content.split('\n'):
        line = line.upper()
        terms = line.split()
        if len(terms) == 0:
            continue

        opcode = terms[0]

        match opcode:
            case "RESIDUE":
                AA_name = terms[1]
                #when we come to ace we stop, tbh idk hat to do here
                if AA_name == "ACE": break
            
            case "ATOM":
                name, type= terms[1], terms[2]
                name_to_type[name] = type[5:]
                if "O" not in name and "S" not in name: seq_set.add(name)
                seq += name
                
            
            case "BOND":
                n = len(terms)
                for i in range(0, n, 3):
                    if terms[i] == "!bond": continue
                    a1, a2 = terms[i+1], terms[i+2]
                    bond = (name_to_type[a1], name_to_type[a2])
                    bond = tuple(sorted(bond))
                    lenght = bond_lenghts[bond]
                    bonds.append((a1, a2, lenght))

            case "END":
                seq_to_AA[seq] = (AA_name, bonds)
                AA_to_seq[AA_name] = seq_set
                bonds = []
                name_to_type = dict()
                seq = ""
                seq_set = set()
                AA_name = ""
            case _:
                pass
    return seq_to_AA, AA_to_seq
           
def compute_dists(atoms : AtomSet, generic_dists : PairAssignment) -> PairAssignment :
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

def write_data(atoms: AtomSet, rhos: TripletAssignment,  cov_dists: PairAssignment, filename = "NOE_data.dat"):

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




            




