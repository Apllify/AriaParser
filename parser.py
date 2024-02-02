ChemShiftToID = dict[tuple[float, float], int]
IDToMol = dict[int, str]

TripletAssignment = list[ tuple[tuple[str, str, str], float] ]

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
        terms = [term.strip for term in line.split(" ")]

        #get all info from current line
        id = int(terms[0])
        chem_shift = float(terms[1])
        err =  float(terms[2])

        mol_type = terms[3]
        res_id = int(terms[4])

        #remember this molecule's name and ID for later
        mol_name =  f"{mol_type}_{res_id}"
        IDToMol[id] = mol_name

        #also store the chem shift of molecule if it is well defined
        if chem_shift != 999 : 
            chem_shift_to_id[(chem_shift - err, chem_shift + err)] = id

        
    return (chem_shift_to_id, id_to_mol)

        





def parse_peaks(content : str, chem_shift_to_id : ChemShiftToID, id_to_mol : IDToMol) -> TripletAssignment:
    pass