"""
Create a ChimeraX defattr file for 293 mix cell entry data
Formatted for the 6ii9_JA_renumbered structure
"""

import pandas as pd
import numpy as np

# Read in the mix functional effects data
func_data = pd.read_csv('results/func_effects/averages/293_mix_entry_func_effects.csv')

# Read in the site numbering map
site_map = pd.read_csv('data/site_numbering_map.csv')

# Calculate average functional effect per site
avg_effects = func_data.groupby('site')['effect'].mean().reset_index()
avg_effects.columns = ['reference_site', 'avg_func_effect']

# Merge with site map to get sequential site numbers and chains
site_data = pd.merge(
    site_map[['reference_site', 'sequential_site', 'chain']],
    avg_effects,
    on='reference_site',
    how='left'
)

# Sort by sequential site
site_data = site_data.sort_values('sequential_site')

# Create the defattr file
output_file = '6ii9_renumbered_mix_func_effects.defattr'

with open(output_file, 'w') as f:
    # Write header
    f.write("# Average functional effect per site\n")
    f.write("# Generated for renumbered 6ii9_JA structure\n")
    f.write("# 293 mix cell entry data\n")
    f.write("# HA2 uses sequential numbering (322-505)\n")
    f.write("# Attribute: avg_func_effect\n")
    f.write("attribute: avg_func_effect\n")
    f.write("match mode: any\n")
    f.write("recipient: residues\n")
    f.write("\n")

    # Write data for each site
    for _, row in site_data.iterrows():
        seq_site = int(row['sequential_site'])
        chain = row['chain']
        effect = row['avg_func_effect']

        if pd.notna(effect):
            f.write(f"\t/{chain}:{seq_site}\t{effect:.4f}\n")

print(f"Defattr file created: {output_file}")
print(f"Total sites with data: {site_data['avg_func_effect'].notna().sum()}")
