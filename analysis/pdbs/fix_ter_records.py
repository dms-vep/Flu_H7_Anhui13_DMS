#!/usr/bin/env python3
"""Fix TER records to reflect actual last residue of each chain."""

# Read the PDB file
with open('6ii9_JA_renumbered.pdb', 'r') as f:
    lines = f.readlines()

# Find the last ATOM record for each chain
last_atom = {}
for line in lines:
    if line.startswith('ATOM'):
        chain = line[21]
        res_num = int(line[22:26].strip())
        res_name = line[17:20].strip()
        atom_serial = int(line[6:11].strip())
        last_atom[chain] = (atom_serial, res_name, res_num)

print("Last atom for each chain:")
for chain, (serial, res_name, res_num) in sorted(last_atom.items()):
    print(f"Chain {chain}: atom {serial}, {res_name} {res_num}")

# Update TER records
new_lines = []
for line in lines:
    if line.startswith('TER'):
        # Get chain from TER record (column 22)
        chain = line[21] if len(line) > 21 else ' '
        if chain in last_atom:
            serial, res_name, res_num = last_atom[chain]
            # TER record format: TER followed by atom serial, residue name, chain, and residue number
            # Columns: 1-6 = "TER   ", 7-11 = serial, 18-20 = resName, 22 = chainID, 23-26 = resSeq
            new_line = f"TER   {serial:5d}      {res_name:>3s} {chain}{res_num:4d}\n"
            new_lines.append(new_line)
            print(f"Updated TER for chain {chain}: {new_line.strip()}")
        else:
            new_lines.append(line)
    else:
        new_lines.append(line)

# Write the updated PDB file
with open('6ii9_JA_renumbered.pdb', 'w') as f:
    f.writelines(new_lines)

print("\nTER records updated successfully!")
