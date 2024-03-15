def extract_residues_in_range(input_file, output_file, x, y):
    """
    Read a prot file and write a new prot file that contains residue in range of [x, y]
    """
    with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
        atom_id = 1  # Start new atom ID from 1
        for line in infile:
            parts = line.split()  # Split line into parts based on whitespace
            if len(parts) >= 5:  # Ensure there are enough columns
                try:
                    residue_number = int(parts[4])  # Fifth column is residue number
                    if x <= residue_number <= y:  # Check if residue number is in the range
                        # Adjust residue number to start from 1
                        new_residue_number = residue_number - x + 1
                        # Write the new line with renumbered atom ID and residue number
                        # Maintaining the spacing between columns similar to the input file
                        new_line = f'{atom_id:>4} {parts[1]:>7} {parts[2]:>5} {parts[3]:<4} {new_residue_number:>4}\n'
                        outfile.write(new_line)
                        atom_id += 1  # Increment new atom ID
                except ValueError:
                    # If conversion to integer fails, skip the line
                    continue

x = 2
y = 3
extract_residues_in_range(f"./data/hmqcnoe.prot", f"./data/hmqcnoe_{x}_{y}.prot", x, y)