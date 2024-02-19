import math
from collections import OrderedDict
from pprint import pprint #debugging tool

#Type defs
Atom = str #  ID of the form ATOMTYPE_RESNUM

ChemShiftToAtom = dict[tuple[float, float], Atom]
AtomsByRes = dict[int, set[Atom]] #atoms by residue
ResIdToAA = dict[int, str]

PairAssignment = dict[ tuple[Atom, Atom], float ]
TripletAssignment = dict[ tuple[Atom, Atom, Atom], float]

PairAngleAssignment = dict[ tuple[Atom, Atom], tuple[float, float, float]] 
TripletAngleAssignment = dict[ tuple[Atom, Atom, Atom], tuple[float, float, float]] 

NOEAssignment = list[ tuple[ list[Atom], list[Atom], list[Atom], float ] ]

ResInfo = tuple[set[Atom], set[Atom],  PairAssignment,  PairAngleAssignment] 
#        (res_content, nonhydrogen_content, res_cov_lengths, res_angle_lengths)
ResInfoDict = dict[str, ResInfo]




#GLOBAL DO NOT TOUCH
CACHE_SIZE = 30
lookup_cache = OrderedDict()


def parse_prot(content : str, res_info_dict : ResInfoDict) -> tuple[ChemShiftToAtom, AtomsByRes, ResIdToAA]:
    """
    Transforms a protein info file in to a dictionary
    which maps each chem shift interval : (LB, UB) to 
    a protein id 

    Ignores all Q atoms for now
    """

    chem_shift_to_atom = dict()
    atom_set = dict()
    res_id_to_AA = dict()

    lines = content.split("\n")
    last_res_id = -1

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

        if last_res_id == -1:
            last_res_id = res_id

        #ignore Qs
        if atom_type[0] == "Q":
            continue

        atom_name =  f"{atom_type.upper()}_{res_id}"
        
        #also store the chem shift of atom if it is well defined
        if chem_shift != 999: 
            chem_shift_to_atom[(chem_shift - err, chem_shift + err)] = atom_name


        #find our AA if we just switch to a new residue
        if last_res_id != res_id:
            AA_name = match_AA(seq, res_info_dict)

            if AA_name == "CYS/SER": 
                res_id_to_AA[last_res_id] = "XAA"
            else: 
                res_id_to_AA[last_res_id] = AA_name

            if AA_name == "MET":
                print(seq)
            atom_set[last_res_id].add(f'O_{last_res_id}')
            #cov_dists[(f'C_{last_res_id}', f'N_{res_id}')] = 1.329 #lenght of peptide bond from aria.par
            
            if AA_name != "XAA": 
                non_hydrogens = res_info_dict[AA_name][1]
                for a in non_hydrogens:
                    seq.add(a)
                    atom_set[last_res_id].add(f'{a}_{last_res_id}')

            last_res_id = res_id
            seq = set()
            hydrogen_counter = 0
        

        #some atom names need to be tweaked to match .top format
        elif atom_type[-1].isdigit() and atom_type[0] == "H" and f'C{atom_type[1:]}' not in seq and f'N{atom_type[1:]}' not in seq:
            hydrogen_counter += 1
            atom_type = f'{atom_type[:-1]}{hydrogen_counter}'
            seq.add(atom_type)
        else:
            hydrogen_counter = 0
            seq.add(atom_type)
        
        atom_name = f'{atom_type}_{res_id}'
        #remember this atom's name and ID for later
        if atom_set.get(res_id) : 
            atom_set[res_id].add(atom_name)
        else : 
            atom_set[res_id] = {atom_name}

        #also store the chem shift of atom if it is well defined
        if chem_shift != 999: 
            chem_shift_to_atom[(chem_shift - err, chem_shift + err)] = atom_name


 
    AA_name = match_AA(seq, res_info_dict)
    if AA_name == "CYS/SER": 
        res_id_to_AA[res_id] = "XAA"        
    else: 
        res_id_to_AA[res_id] = AA_name

    if AA_name != "XAA": 
        non_hydrogens = res_info_dict[AA_name][1]
        for a in non_hydrogens:
            atom_set[res_id].add(f'{a}_{last_res_id}')

            
    return (chem_shift_to_atom, atom_set, res_id_to_AA)

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
    Tries to find the amino acid which corresponding to resiude.
    
    Returns "XAA" if no match found
    """
    for name, (seq,_, _, _) in res_info_dict.items():
        if seq == atoms:
           #special case since CYS and SER share same composition when removing hydrogens
            if name in ("CYS", "SER"):
               return "CYS/SER" 
           
            return name
    
    if atoms | {"HZ3"} == res_info_dict["LYS"][0]: return "LYS"
    if atoms | {"HE2"} == res_info_dict["HIS"][0]: return "HIS"
    if atoms - {"HE1"} == res_info_dict["GLU"][0]: return "GLU"
    if atoms - {"HD1"} == res_info_dict["ASP"][0]: return "ASP"
    if {"NE", "NH1", "NH2"} <= atoms: return "ARG"

    return "XAA"
    
def parse_peaks(content : str, chem_shift_to_atom : ChemShiftToAtom) -> NOEAssignment:
    """
    Uses the prot information to read the peaks file 
    and transform each triplet of chem shifts into their 
    respective atom IDs and stores the NOE measure for each.
    """
    
    lines = content.split("\n")
    #Assignmnet of rho to corresponding atoms follows scheme detailed in README.txt
    triplet_assignment : NOEAssignment = []


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
        Iq = list(filter(lambda atom : atom[0] == "H", Iq))
        Ir = list(filter(lambda atom : atom[0] == "H", Ir))
    
        #Check that each chemical shift correspnds to something
        if 0 in (len(Ip), len(Iq), len(Ir)):
            continue
    

        triplet_assignment.append((Ip, Iq, Ir, noe))

    return triplet_assignment

def parse_par(content : str) -> tuple[PairAssignment, TripletAngleAssignment]: 
    """
    Parses the content of the .par file
    into a list of generic distances between atom types

    Returns a tuple where elem 1 is the assignment obtained from 
    bond statements, while elem2 is the assignment obtained
    from angle statements.
    """

    bond_assignment : PairAssignment = dict()
    angle_assignment : TripletAngleAssignment = dict()

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
                pair_key = tuple(sorted((m1, m2)))
                bond_assignment[pair_key] = dist

            case "ANGLE" : 
                #retrieve angle data
                try : 
                    m1, m2, m3 = terms[1], terms[2], terms[3]
                    angle = float(terms[5])
                except : 
                    continue

                #get the two lengths of the triangle
                edge1, edge2 = tuple(sorted((m1, m2))), tuple(sorted((m2, m3)))
                d1 = bond_assignment.get(edge1)
                d2 = bond_assignment.get(edge2)
                if d1 == None or d2 == None : 
                    continue

                #compute and store length + extra info on this triplet
                dist = math.sqrt( d1**2 + d2**2 - d1*d2*math.cos(math.radians(angle)) ) 
                d1d2 = d1*d2
                sin_alpha = math.sin(math.radians(angle))

                pair_key = (m2, *sorted((m1, m3))) #anchor atom first

                angle_assignment[pair_key] = (dist, d1d2, sin_alpha) #need a lot of info for model slack

    return (bond_assignment, angle_assignment)

def parse_top(content: str, cov_lengths: PairAssignment, angle_lengths : TripletAngleAssignment) -> ResInfoDict :
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
    res_cov_lengths : PairAssignment = dict()
    res_angle_lengths : PairAngleAssignment = dict()

    res_atom_to_type = dict() #doesn't appear in output
    non_hydrogens = set()


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

                    #keep two diff running lists of atoms
                    if name[0] not in ("O", "S"): 
                        res_atoms.add(name)
                    else:
                        non_hydrogens.add(name)
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

                        #check if known cov bond length
                        cov_length = cov_lengths.get(general_bond)
                        if cov_length != None :
                            res_cov_lengths[(a1, a2)] = cov_length
                            continue

                case "END":

                    #go through atom triplets to check if any dists are known from angles
                    res_bonds = list(res_cov_lengths.keys())
                    for anchor in res_atoms:  
                        for i in range(len(res_bonds)):

                            #consider a first bond of our central atom
                            bond_i = res_bonds[i]
                            if anchor not in bond_i : 
                                continue 
                            try : 
                                atom_i = next(filter(lambda x : x != anchor, bond_i))
                            except : 
                                continue

                            for j in range(i+1, len(res_bonds)):
                                #consider a second different bond to our central atom
                                bond_j = res_bonds[j]
                                if anchor not in bond_j : 
                                    continue 
                                try : 
                                    atom_j = next(filter(lambda x : x != anchor, bond_j))
                                except : 
                                    continue


                                #check if this triplet is in the angles assignment
                                general_ij = (res_atom_to_type[atom_i], res_atom_to_type[atom_j])
                                triplet_key = (res_atom_to_type[anchor], *sorted(general_ij))


                                #store this length if applicable
                                angle_info  =angle_lengths.get(triplet_key)
                                if angle_info != None :

                                    res_angle_lengths[(atom_i, atom_j)] = angle_info

                    #add our new entry and reset variables
                    res_info_dict[res_name] = (res_atoms, non_hydrogens, res_cov_lengths, res_angle_lengths)

                    inside_res = False
                    res_name = ""
                    res_atoms = set()
                    non_hydrogens = set()

                    res_cov_lengths = dict()
                    res_angle_lengths = dict()

                    res_atom_to_type = dict()
        
        #merge CYS and SER information
        if res_info_dict.get("SER") and res_info_dict.get("CYS"):
            a_ser, _, cov_ser, ang_ser = res_info_dict["SER"]
            a_cys, _, cov_cys, ang_cys = res_info_dict["SER"]

            res_info_dict["CYS/SER"] = (a_ser & a_cys, [], dict(cov_ser.items() & cov_cys.items()),  dict(ang_ser.items() & ang_cys.items()))
    
    return res_info_dict
           

def compute_dists(atoms : AtomsByRes, res_info : ResInfoDict, res_id_to_AA : ResIdToAA) -> tuple[PairAssignment, PairAssignment] :
    """
    Uses the generic bond + angle distances 
    obtained from the .par + .top file to find all the
    distances applicable to our atoms set

    Returns a tuple of the form (covalent_dists, angle_dists)
    """


    cov_dists : PairAssignment = dict()
    ang_dists : PairAssignment = dict()
    last_res_id = -10

    #proceed by residue
    for res_id, residue in atoms.items():
        
        cur_AA = res_id_to_AA[res_id]
        
        #ignore unkonwn residue types
        if cur_AA == "XAA":
            continue

        #hard-coded length of peptide bond
        if (res_id - last_res_id) == 1: 
            cov_dists[(f'C_{last_res_id}', f'N_{res_id}')] = 1.329 

        #go through all angle + bond lengths
        _, _, cov_lengths, ang_lengths = res_info[cur_AA]
        for (a1, a2), length in cov_lengths.items()  :
            a1_spec = f"{a1}_{res_id}"
            a2_spec = f"{a2}_{res_id}"

            if {a1_spec, a2_spec} <= residue : 
                cov_dists[(a1_spec, a2_spec)] = length

        for (a1, a2), length in ang_lengths.items()    :
            a1_spec = f"{a1}_{res_id}"
            a2_spec = f"{a2}_{res_id}"

            if {a1_spec, a2_spec} <= residue : 
                ang_dists[(a1_spec, a2_spec)] = length

        last_res_id = res_id
            
    return (cov_dists, ang_dists)



def write_data(atoms: AtomsByRes, rhos: NOEAssignment,  cov_dists: PairAssignment, ang_dists : PairAssignment, filename = "NOE_data.dat"):

    with open(filename, "w") as outfile:
        #define atoms set
        outfile.write("set ATOMS := ")
        residues = atoms.values()
        for residue in residues : 
            for a in residue: outfile.write(f' {a}')
        outfile.write(";\n")

        #define angle distances set 
        outfile.write("set ANGDISTS :=")
        for a1, a2 in ang_dists.keys() : outfile.write(f" {a1} {a2}")
        outfile.write(";\n")

        #define covalent distances set
        outfile.write("set COVDISTS := ")
        for a1, a2 in cov_dists.keys(): outfile.write(f' {a1} {a2}')
        outfile.write(";\n")

        #define RHOS set
        num_rhos = len(rhos)
        outfile.write(f'param numRHOS := {num_rhos};\n')

        #define NOE_A1 and A_2 sets:
        i = 1
        for _, Iq, Ir, rho in rhos:
            outfile.write(f'set NOE_A1[{i}] := ')
            for a in Iq:
                outfile.write(f'{a} ')
            outfile.write(";\n")
            outfile.write(f'set NOE_A2[{i}] := ')
            for a in Ir:
                outfile.write(f'{a} ')
            outfile.write(";\n")
            i += 1

        #give covalent distance data
        outfile.write("param  CovDists := ")
        for (a1, a2), dist in cov_dists.items(): outfile.write(f' {a1} {a2} {dist}')
        outfile.write(";\n")

        #give angle distance data
        outfile.write("param  AngDists := ")
        for (a1, a2), (dist, _, _) in ang_dists.items(): outfile.write(f' {a1} {a2} {dist}')
        outfile.write(";\n")
        outfile.write("param  d1d2 := ")
        for (a1, a2), (_, d1d2, _) in ang_dists.items(): outfile.write(f' {a1} {a2} {d1d2}')
        outfile.write(";\n")
        outfile.write("param  SinAlpha := ")
        for (a1, a2), (_, _, sin_alpha) in ang_dists.items(): outfile.write(f' {a1} {a2} {sin_alpha}')
        outfile.write(";\n")

        #give RHO data
        i = 1
        outfile.write("param rho := ")
        for _, _, _, rho in rhos:
            outfile.write(f'{i} {rho} ')
            i += 1
        outfile.write(";\n")