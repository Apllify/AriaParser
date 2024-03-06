import utils.atom_class as atom_class

def get_correct_atom_name(residue, res_id, atoms_to_coord):
    """
    Convert atom names from H2 H3 to H1 H2, and H to HN
    """
    residue_atoms = [atom_class.Atom(atom.name, atom.coord) for atom in residue]
    residue_atoms = sorted(residue_atoms)
    for i in range(len(residue_atoms)):
        # Change all 2-3 hydrogen atoms to 1-2 hydrogen atoms
        if i < len(residue_atoms)-2 \
        and residue_atoms[i+1].name[0] == 'H' and residue_atoms[i+1].name[-1] == '2' \
        and residue_atoms[i+2].name[0] == 'H' and residue_atoms[i+2].name[-1] == '3' \
        and residue_atoms[i].name != residue_atoms[i+1].name[:-1] + '1':
            residue_atoms[i+1].name = residue_atoms[i+1].name[:-1] + '1'
            residue_atoms[i+2].name = residue_atoms[i+2].name[:-1] + '2'
        # Change all H atoms to HN atoms
        if residue_atoms[i].name == 'H':
            residue_atoms[i].name = 'HN'
        # Add res_id to name
        residue_atoms[i].name = f'{residue_atoms[i].name}_{res_id}'
    return sorted([atom for atom in residue_atoms if atom.name in atoms_to_coord])

def get_aa_names_from_file(file):
    """
    Get all aa names in a top file
    """
    with open(f"{file}", "r") as stream:
        content = stream.read()
    aa_names = []
    res_name = ""

    for line in content.split('\n'):
        line = line.upper()
        terms = line.split()
        if len(terms) == 0:
            continue
        opcode = terms[0]
        if opcode == "RESIDUE": 
            res_name = terms[1]
            if res_name == "ACE": 
                break #rest of file is irrelevant 
            if res_name not in ['ZN', 'CHEX', 'HISH']:
                aa_names.append(res_name)
    return aa_names