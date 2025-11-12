#!/usr/bin/env python3
"""
Convert functional effects CSV to ChimeraX defattr format for 6ii9 structure.
This structure has HA2 numbering that restarts at 1.
"""

import pandas as pd
import sys

# Configuration
input_file = "scratch_notebooks/dmsviz/2-6_entry_func_effects_sequential.csv"
output_file = "6ii9_2-6_func_effects.defattr"
# For 6ii9: HA2 numbering restarts, so we need to adjust
# HA1 ends at sequential site 317 (cleavage site at 318-321)
# HA2 starts at sequential site 323
ha2_start_sequential = 323  # First HA2 residue in sequential numbering
ha2_restart_at = 1  # HA2 numbering starts at 1 in 6ii9

print(f"Reading {input_file}...")
df = pd.read_csv(input_file)

print(f"Data contains {len(df)} rows")
print(f"Sites: {df['site'].min()} to {df['site'].max()}")

# Calculate average effect per site (across all mutations)
site_avg_effect = df.groupby('site')['effect'].mean().reset_index()
site_avg_effect.columns = ['site', 'avg_effect']

# Sort by site to ensure correct order
site_avg_effect = site_avg_effect.sort_values('site').reset_index(drop=True)

print(f"\nCalculated average effects for {len(site_avg_effect)} sites")
print(f"Sequential sites: {site_avg_effect['site'].min()} to {site_avg_effect['site'].max()}")

# Create PDB numbering for 6ii9
# HA1 region: sequential sites 1-317 map directly
# Cleavage site: sequential sites 318-321 (skip these for now)
# HA2 region: sequential sites 323+ map to PDB numbering starting at 1
def get_pdb_residue_number(sequential_site):
    if sequential_site < ha2_start_sequential:
        # HA1 region - use sequential numbering directly
        return int(sequential_site)
    else:
        # HA2 region - renumber starting from 1
        return int(sequential_site - ha2_start_sequential + ha2_restart_at)

# Add PDB numbering column
site_avg_effect['pdb_residue'] = site_avg_effect['site'].apply(get_pdb_residue_number)

print(f"\nHA2 starts at sequential site {ha2_start_sequential}")
print(f"HA2 will be numbered starting at {ha2_restart_at} in the PDB")
print(f"Example: Sequential site {ha2_start_sequential} -> PDB residue {get_pdb_residue_number(ha2_start_sequential)}")

# Write defattr file
# For trimer: typically chains are (A,B), (C,D), (E,F) for HA1,HA2 pairs
chain_pairs = [('A', 'B'), ('C', 'D'), ('E', 'F')]

with open(output_file, 'w') as f:
    # Write header
    f.write("# Average functional effect per site\n")
    f.write("# Generated for 6ii9 structure (HA2 numbering restarts at 1)\n")
    f.write("# Attribute: avg_func_effect\n")
    f.write("attribute: avg_func_effect\n")
    f.write("match mode: any\n")
    f.write("recipient: residues\n")
    f.write("\n")

    # Write data lines for each protomer
    for ha1_chain, ha2_chain in chain_pairs:
        for _, row in site_avg_effect.iterrows():
            sequential_site = row['site']
            pdb_residue = row['pdb_residue']
            effect = row['avg_effect']

            # Determine chain based on region
            if sequential_site < ha2_start_sequential:
                # HA1 region
                chain = ha1_chain
            else:
                # HA2 region
                chain = ha2_chain

            f.write(f"\t/{chain}:{int(pdb_residue)}\t{effect:.4f}\n")

print(f"\nWrote {output_file}")
print(f"\nTo use in ChimeraX:")
print(f"  1. Open your 6ii9 structure")
print(f"  2. Run: open {output_file}")
print(f"  3. Color by attribute: color byattribute avg_func_effect palette red:white")
print(f"\nNote: This file includes all 3 protomers with chains:")
print(f"  Protomer 1: A (HA1), B (HA2)")
print(f"  Protomer 2: C (HA1), D (HA2)")
print(f"  Protomer 3: E (HA1), F (HA2)")
