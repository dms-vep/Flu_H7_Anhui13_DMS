"""
Create a ChimeraX defattr file for 293 mix cell entry data
Formatted for the 6ii9_JA_renumbered structure
Using FILTERED data (times_seen >= 2)

Maps from reference_site numbering (H3 numbering) to sequential_site numbering (PDB numbering)
"""

import pandas as pd
import numpy as np
import re

# Read in the FILTERED cell entry data
# This data has already been filtered for times_seen >= 2
# Site numbers are in reference_site numbering (H3 numbering)
func_data = pd.read_csv('filtered_cell_entry.csv')

# Calculate average functional effect per site from FILTERED data
avg_effects = func_data.groupby('site')['effect'].mean().reset_index()
avg_effects.columns = ['reference_site', 'avg_func_effect']

# Read the site numbering map to convert reference_site to sequential_site (PDB numbering)
site_map = pd.read_csv('../../../data/site_numbering_map.csv')

# Both filtered_cell_entry and site_map use reference_site as strings
# Keep them as strings for merging
avg_effects['reference_site_str'] = avg_effects['reference_site'].astype(str)

# Also ensure site_map reference_site is string
site_map['reference_site_str'] = site_map['reference_site'].astype(str)

# Merge with site map on reference_site (as strings)
site_data = pd.merge(
    avg_effects,
    site_map[['sequential_site', 'reference_site_str', 'chain']],
    left_on='reference_site_str',
    right_on='reference_site_str',
    how='left'
)

# For insertion sites (like 158a, 158b, 276a, 325a), they won't match in the merge
# because the site_map only has non-insertion sites
sites_with_mapping = site_data[site_data['sequential_site'].notna()].copy()
sites_without_mapping = site_data[site_data['sequential_site'].isna()].copy()

print(f"Sites with sequential mapping: {len(sites_with_mapping)}")
if len(sites_without_mapping) > 0:
    print(f"Sites without sequential mapping (likely insertions): {len(sites_without_mapping)}")
    print(f"  These sites: {sites_without_mapping['reference_site_str'].tolist()[:20]}...")  # Show first 20

# Calculate statistics for the header
min_effect = sites_with_mapping['avg_func_effect'].min()
max_effect = sites_with_mapping['avg_func_effect'].max()

# Read the PDB file to get actual residue ranges per chain
# This ensures we only include residues that exist in the structure
pdb_file = '../../../analysis/pdbs/6ii9_JA_renumbered.pdb'
pdb_residues = {}
with open(pdb_file, 'r') as f:
    for line in f:
        if line.startswith('ATOM'):
            chain = line[21]
            resnum = int(line[22:27].strip())
            if chain not in pdb_residues:
                pdb_residues[chain] = set()
            pdb_residues[chain].add(resnum)

# Get ranges for each chain
for chain in pdb_residues:
    pdb_residues[chain] = sorted(pdb_residues[chain])

print(f"\nPDB structure residue ranges:")
for chain in sorted(pdb_residues.keys()):
    print(f"  Chain {chain}: {min(pdb_residues[chain])}-{max(pdb_residues[chain])} ({len(pdb_residues[chain])} residues)")

# Separate HA1 and HA2 based on sequential_site
# HA1 chains A,C,E: filter by actual PDB range
# HA2 chains B,D,F: filter by actual PDB range
ha1_max = max(pdb_residues.get('A', [0]) + pdb_residues.get('C', [0]) + pdb_residues.get('E', [0]))
ha2_min = min(pdb_residues.get('B', [999]) + pdb_residues.get('D', [999]) + pdb_residues.get('F', [999]))
ha2_max = max(pdb_residues.get('B', [0]) + pdb_residues.get('D', [0]) + pdb_residues.get('F', [0]))

ha1_data = sites_with_mapping[sites_with_mapping['sequential_site'] <= ha1_max].copy()
ha2_data = sites_with_mapping[(sites_with_mapping['sequential_site'] >= ha2_min) &
                               (sites_with_mapping['sequential_site'] <= ha2_max)].copy()

print(f"\nFiltering to PDB ranges:")
print(f"  HA1: 1-{ha1_max}")
print(f"  HA2: {ha2_min}-{ha2_max}")

# Create the defattr file
output_file = '../../../chimera_293_mix_entry.defattr'

with open(output_file, 'w') as f:
    # Write header with metadata
    f.write("# ChimeraX attribute file for H7 Anhui13 cell entry effects in 293 mix cells\n")
    f.write("# Data: Site-level average effects (FILTERED: times_seen >= 2)\n")
    f.write(f"# Effect range: {min_effect:.3f} to {max_effect:.3f}\n")
    f.write("# Color scheme: redblue (Vega diverging red-blue scheme)\n")
    f.write("#   Red (#67001f) = negative effects (reduced entry)\n")
    f.write("#   Blue (#053061) = positive effects (enhanced entry)\n")
    f.write("#\n")
    f.write("# Filtering applied: times_seen >= 2\n")
    f.write("# Residue numbering: SEQUENTIAL (matches 6ii9_JA_renumbered.pdb)\n")
    f.write("#   HA1 chains A,C,E: sequential sites 1-317\n")
    f.write("#   HA2 chains B,D,F: sequential sites 322-491\n")
    f.write("# Applied to all protomers in the trimer structure\n")
    f.write("#\n")
    f.write("# To apply colors in ChimeraX, use:\n")
    f.write(f"#   color byattribute cell_entry_effect palette #67001f:#053061 range {min_effect:.3f},{max_effect:.3f}\n")
    f.write("\n")
    f.write("attribute: cell_entry_effect\n")
    f.write("recipient: residues\n")
    f.write("match mode: any\n")
    f.write("\n")

    # Write HA1 data (chains A, C, E)
    # Use sequential_site numbering for the PDB
    # Check each residue exists in each chain before writing
    for _, row in ha1_data.iterrows():
        seq_site = int(row['sequential_site'])
        effect = row['avg_func_effect']

        if pd.notna(effect):
            # Apply to HA1 protomers, but only if residue exists in that chain
            for chain in ['A', 'C', 'E']:
                if seq_site in pdb_residues.get(chain, []):
                    f.write(f"\t/{chain}:{seq_site}\t{effect:.4f}\n")

    # Write HA2 data (chains B, D, F)
    # Use sequential_site numbering for the PDB
    # Check each residue exists in each chain before writing
    for _, row in ha2_data.iterrows():
        seq_site = int(row['sequential_site'])
        effect = row['avg_func_effect']

        if pd.notna(effect):
            # Apply to HA2 protomers, but only if residue exists in that chain
            for chain in ['B', 'D', 'F']:
                if seq_site in pdb_residues.get(chain, []):
                    f.write(f"\t/{chain}:{seq_site}\t{effect:.4f}\n")

print(f"Defattr file created: {output_file}")
print(f"Total HA1 sites with data: {len(ha1_data)}")
print(f"Total HA2 sites with data: {len(ha2_data)}")
print(f"Total sites: {len(sites_with_mapping)}")
print(f"Effect range: {min_effect:.4f} to {max_effect:.4f}")
if len(sites_without_mapping) > 0:
    print(f"\nWarning: {len(sites_without_mapping)} sites could not be mapped")
    unmapped_sites = sorted(sites_without_mapping['reference_site_str'].unique().tolist())
    print(f"  These reference sites will be skipped: {unmapped_sites[:50]}...")  # Show first 50
    # Check which are likely insertions
    insertions = [s for s in unmapped_sites if not s.isdigit()]
    if insertions:
        print(f"  Insertions (not in PDB): {insertions}")
