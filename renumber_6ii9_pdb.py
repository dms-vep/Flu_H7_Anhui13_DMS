#!/usr/bin/env python3
"""
Renumber 6ii9_JA.pdb file to match sequential numbering scheme.
- HA1 stays the same (residues 1-317)
- HA2 currently numbered 1-184 in PDB, will be renumbered to 322-505
  (Sequential 323-506, minus 1 because cleavage site PKGR is missing)
"""

import sys

input_pdb = "scratch_notebooks/dmsviz/6ii9_JA.pdb"
output_pdb = "6ii9_JA_renumbered.pdb"

# Configuration
# HA2 chains in the PDB (assuming B, D, F based on typical trimer structure)
ha2_chains = ['B', 'D', 'F']
# HA2 starts at 1 in original PDB, we want it to start at 322
# Sequential site 323 (first HA2 residue 'G') will be 322 in renumbered PDB
# because the 4 aa cleavage site (318-321) is missing
ha2_offset = 321  # Add this to original HA2 residue numbers

print(f"Reading {input_pdb}...")
print(f"HA2 chains: {ha2_chains}")
print(f"HA2 residues will be offset by +{ha2_offset}")
print(f"Example: Original /B:1 -> Renumbered /B:322")

with open(input_pdb, 'r') as infile, open(output_pdb, 'w') as outfile:
    for line in infile:
        # Process ATOM and HETATM records
        if line.startswith(('ATOM', 'HETATM')):
            # PDB format: residue number is in columns 23-26 (0-indexed: 22-26)
            chain = line[21]

            if chain in ha2_chains:
                # Extract original residue number
                try:
                    old_resnum = int(line[22:26].strip())
                    new_resnum = old_resnum + ha2_offset

                    # Replace residue number in the line
                    # Format as right-justified 4-character field
                    new_line = line[:22] + f"{new_resnum:4d}" + line[26:]
                    outfile.write(new_line)
                except ValueError:
                    # If residue number has insertion code or other issue, keep original
                    outfile.write(line)
            else:
                # HA1 or other chains - keep original
                outfile.write(line)
        else:
            # Non-ATOM lines (headers, etc.) - keep as is
            outfile.write(line)

print(f"\nWrote renumbered PDB to: {output_pdb}")
print(f"\nTo use with ChimeraX:")
print(f"  1. Open the renumbered PDB: open {output_pdb}")
print(f"  2. Use the original defattr file (293_2-6_func_effects_sequential.defattr)")
print(f"     which uses sequential numbering")
