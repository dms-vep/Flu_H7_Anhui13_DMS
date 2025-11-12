#!/usr/bin/env python3
"""
Convert functional effects CSV to ChimeraX defattr format for renumbered 6ii9 structure.
This version uses sequential numbering throughout (HA2 continues from 322, not restart at 1).
"""

import pandas as pd

# Configuration
input_file = "scratch_notebooks/dmsviz/2-6_entry_func_effects_sequential.csv"
output_file = "6ii9_renumbered_2-6_func_effects.defattr"

# For renumbered 6ii9_JA.pdb:
# HA1: chains A, C, E (residues 1-317)
# HA2: chains B, D, F (residues 322-505, sequential 323-506 minus cleavage site)
# The cleavage site (sequential 318-321) is missing in the PDB

# Chain mapping for the 3 protomers
chain_pairs = [('A', 'B'), ('C', 'D'), ('E', 'F')]

print(f"Reading {input_file}...")
df = pd.read_csv(input_file)

print(f"Data contains {len(df)} rows")
print(f"Sites: {df['site'].min()} to {df['site'].max()}")

# Calculate average effect per site (across all mutations)
site_avg_effect = df.groupby('site')['effect'].mean().reset_index()
site_avg_effect.columns = ['site', 'avg_effect']

# Sort by site
site_avg_effect = site_avg_effect.sort_values('site').reset_index(drop=True)

print(f"\nCalculated average effects for {len(site_avg_effect)} sites")
print(f"Sequential sites: {site_avg_effect['site'].min()} to {site_avg_effect['site'].max()}")

# Map sequential sites to PDB residue numbers for renumbered PDB
def get_pdb_residue_number(sequential_site):
    if sequential_site <= 317:
        # HA1 region (1-317)
        return int(sequential_site)
    elif 318 <= sequential_site <= 321:
        # Cleavage site - not in PDB, skip
        return None
    else:
        # HA2 region (sequential 323-506 -> PDB 322-505)
        # Sequential 323 -> PDB 322 (offset by -1 because cleavage missing)
        return int(sequential_site - 1)

print("\nMapping sequential to PDB numbering:")
print("  Sequential 1-317 -> PDB 1-317 (HA1)")
print("  Sequential 318-321 -> SKIPPED (cleavage site not in PDB)")
print("  Sequential 323-506 -> PDB 322-505 (HA2, offset by -1)")

# Write defattr file
with open(output_file, 'w') as f:
    # Write header
    f.write("# Average functional effect per site\n")
    f.write("# Generated for renumbered 6ii9_JA structure\n")
    f.write("# HA2 uses sequential numbering (322-505)\n")
    f.write("# Attribute: avg_func_effect\n")
    f.write("attribute: avg_func_effect\n")
    f.write("match mode: any\n")
    f.write("recipient: residues\n")
    f.write("\n")

    # Write data lines for each protomer
    for ha1_chain, ha2_chain in chain_pairs:
        for _, row in site_avg_effect.iterrows():
            sequential_site = row['site']
            pdb_residue = get_pdb_residue_number(sequential_site)

            # Skip cleavage site residues (not in PDB)
            if pdb_residue is None:
                continue

            effect = row['avg_effect']

            # Determine chain based on region
            if sequential_site <= 317:
                # HA1 region
                chain = ha1_chain
            else:
                # HA2 region (sequential >= 323)
                chain = ha2_chain

            f.write(f"\t/{chain}:{int(pdb_residue)}\t{effect:.4f}\n")

print(f"\nWrote {output_file}")
print(f"\nTo use in ChimeraX:")
print(f"  1. Open renumbered structure: open 6ii9_JA_renumbered.pdb")
print(f"  2. Load attributes: open {output_file}")
print(f"  3. Color by attribute: color byattribute avg_func_effect palette red:white")
print(f"\nNote: This file includes all 3 protomers with chains:")
print(f"  Protomer 1: A (HA1), B (HA2)")
print(f"  Protomer 2: C (HA1), D (HA2)")
print(f"  Protomer 3: E (HA1), F (HA2)")
