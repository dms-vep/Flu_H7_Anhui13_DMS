#!/usr/bin/env python3
"""
Convert functional effects CSV to ChimeraX defattr format.
"""

import pandas as pd
import sys

# Configuration
input_file = "scratch_notebooks/dmsviz/2-6_entry_func_effects_sequential.csv"
output_file = "293_2-6_func_effects_sequential.defattr"
# For trimers with repeated chains, specify all chains to apply attributes to
# Set to None to apply without chain specification, or list specific chains
chain_ids = ["A", "B"]  # Will create entries for both chains A and B (covering all protomers)

print(f"Reading {input_file}...")
df = pd.read_csv(input_file)

print(f"Data contains {len(df)} rows")
print(f"Sites: {df['site'].min()} to {df['site'].max()}")

# Calculate average effect per site (across all mutations)
# Site column already contains sequential numbering
site_avg_effect = df.groupby('site')['effect'].mean().reset_index()
site_avg_effect.columns = ['site', 'avg_effect']

# Sort by site to ensure correct order
site_avg_effect = site_avg_effect.sort_values('site').reset_index(drop=True)

print(f"\nCalculated average effects for {len(site_avg_effect)} sites")
print(f"Sequential sites: {site_avg_effect['site'].min()} to {site_avg_effect['site'].max()}")

# Write defattr file
with open(output_file, 'w') as f:
    # Write header
    f.write("# Average functional effect per site\n")
    f.write("# Generated from 293_2-6_entry_func_effects.csv\n")
    f.write("# Attribute: avg_func_effect\n")
    f.write("attribute: avg_func_effect\n")
    f.write("match mode: any\n")  # Use 'any' mode to allow matching multiple residues
    f.write("recipient: residues\n")
    f.write("\n")

    # Write data lines (format: /chain:residue_number    value)
    # If multiple chains specified, write entries for each chain
    for _, row in site_avg_effect.iterrows():
        site = row['site']
        # Convert to int if it's a whole number, otherwise keep as string (for insertion codes)
        if isinstance(site, float) and site.is_integer():
            site = int(site)
        site = str(site)
        effect = row['avg_effect']
        if chain_ids:
            for chain in chain_ids:
                f.write(f"\t/{chain}:{site}\t{effect:.4f}\n")
        else:
            f.write(f"\t:{site}\t{effect:.4f}\n")

print(f"\nWrote {output_file}")
if chain_ids:
    print(f"Attributes assigned to chains: {', '.join(chain_ids)}")
print(f"\nTo use in ChimeraX:")
print(f"  1. Open your structure")
print(f"  2. Run: open {output_file}")
print(f"  3. Color by attribute: color byattribute avg_func_effect palette <palette_name>")
if chain_ids:
    print(f"\nNote: Attributes are assigned to all specified chains.")
    print(f"This will cover all protomers in your trimer structure.")
