#!/usr/bin/env python3
"""
Create a table showing average functional effects with both sequential and reference numbering.
This helps verify that the defattr file is correctly mapping data to the structure.
"""

import pandas as pd

# Read the functional effects data
print("Reading functional effects data...")
func_data = pd.read_csv('scratch_notebooks/dmsviz/2-6_entry_func_effects_sequential.csv')

# Calculate average effect per site
site_avg = func_data.groupby('site')['effect'].mean().reset_index()
site_avg.columns = ['sequential_site', 'avg_func_effect']

# Read site numbering map
print("Reading site numbering map...")
site_map = pd.read_csv('data/site_numbering_map.csv')

# Merge to get both numbering schemes
print("Merging data...")
table = pd.merge(
    site_avg,
    site_map[['sequential_site', 'reference_site', 'sequential_wt', 'region']],
    on='sequential_site',
    how='left'
)

# Calculate PDB numbering for renumbered 6ii9 structure
def get_6ii9_pdb_number(seq_site):
    if seq_site <= 317:
        return seq_site
    elif 318 <= seq_site <= 321:
        return 'N/A (cleavage)'
    else:
        return seq_site - 1

table['6ii9_renumbered_pdb'] = table['sequential_site'].apply(get_6ii9_pdb_number)

# Reorder columns
table = table[['sequential_site', 'reference_site', '6ii9_renumbered_pdb',
               'sequential_wt', 'region', 'avg_func_effect']]

# Sort by sequential site
table = table.sort_values('sequential_site').reset_index(drop=True)

# Save to CSV
output_file = 'site_numbering_verification_table.csv'
table.to_csv(output_file, index=False)

print(f"\nCreated table with {len(table)} sites")
print(f"Saved to: {output_file}")

# Display first few rows
print("\nFirst 20 rows:")
print(table.head(20).to_string(index=False))

print("\n\nHA1/HA2 transition (rows around cleavage site):")
print(table[(table['sequential_site'] >= 315) & (table['sequential_site'] <= 325)].to_string(index=False))

print("\n\nLast 10 rows:")
print(table.tail(10).to_string(index=False))

# Show some statistics
print("\n\nSummary statistics:")
print(f"  Sequential sites: {table['sequential_site'].min()} to {table['sequential_site'].max()}")
print(f"  Avg func effect range: {table['avg_func_effect'].min():.4f} to {table['avg_func_effect'].max():.4f}")
print(f"  Mean effect: {table['avg_func_effect'].mean():.4f}")
print(f"  Sites with strong negative effects (<-4): {(table['avg_func_effect'] < -4).sum()}")
print(f"  Sites with positive effects (>0): {(table['avg_func_effect'] > 0).sum()}")
