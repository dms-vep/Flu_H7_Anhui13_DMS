"""
Create a full-page heatmap of H7 cell entry effects on 2,6 cells with region annotations.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors
import seaborn as sns

# Read in the data
func_data = pd.read_csv('results/func_effects/averages/293_2-6_entry_func_effects.csv')
site_map = pd.read_csv('data/site_numbering_map.csv')

# Merge and filter data
func_data_for_heatmap = pd.merge(
    func_data,
    site_map[['reference_site', 'region', 'rbs_region']],
    left_on='site',
    right_on='reference_site',
    how='right'
).drop(
    columns=['reference_site']
).query(
    'times_seen >= 2'
).query(
    'effect_std <= 2'
)

# Filter to only include amino acid mutations (not stop codons or wildtype)
heatmap_data = func_data_for_heatmap[
    (func_data_for_heatmap['mutant'] != '*') &
    (func_data_for_heatmap['mutant'] != func_data_for_heatmap['wildtype'])
].copy()

# Pivot the data to create a matrix for the heatmap
# Rows = amino acids, Columns = sites
amino_acids = list('ACDEFGHIKLMNPQRSTVWY')
heatmap_matrix = heatmap_data.pivot_table(
    index='mutant',
    columns='site',
    values='effect',
    aggfunc='first'
)

# Reindex to ensure all amino acids are present and in order
heatmap_matrix = heatmap_matrix.reindex(amino_acids)

# Sort columns by site number
heatmap_matrix = heatmap_matrix.sort_index(axis=1)

# Create region annotations for the top
site_regions = site_map.set_index('reference_site')['region'].to_dict()
regions_ordered = [site_regions.get(site, '') for site in heatmap_matrix.columns]

# Define region colors
region_colors = {
    'epitope-A': '#E8B44F',
    'epitope-B': '#E8B44F',
    'epitope-C': '#E8B44F',
    'epitope-D': '#E8B44F',
    'epitope-E': '#E8B44F',
    'HA1': '#D3D3D3',
    'HA2': '#A0A0A0',
    'cleavage site': '#1A391A',
}

# Create the figure
fig = plt.figure(figsize=(24, 8))
gs = fig.add_gridspec(2, 1, height_ratios=[0.5, 20], hspace=0.02)

# Top panel: Region annotations
ax_regions = fig.add_subplot(gs[0])
# Convert hex colors to RGB for imshow
from matplotlib.colors import to_rgb
region_bar = np.array([[to_rgb(region_colors.get(r, '#FFFFFF')) for r in regions_ordered]])
ax_regions.imshow(region_bar, aspect='auto', interpolation='nearest')
ax_regions.set_xticks([])
ax_regions.set_yticks([])
ax_regions.set_ylabel('Region', fontsize=10, rotation=0, ha='right', va='center')

# Add region labels at the top
current_region = None
region_start = 0
for i, region in enumerate(regions_ordered):
    if region != current_region:
        if current_region is not None:
            # Add label for previous region
            mid_point = (region_start + i - 1) / 2
            if 'epitope' in current_region:
                label = current_region.replace('epitope-', '')
                ax_regions.text(mid_point, -0.5, label, ha='center', va='bottom',
                               fontsize=12, fontweight='bold')
        current_region = region
        region_start = i

# Add last region label
if current_region is not None and region_start < len(regions_ordered):
    mid_point = (region_start + len(regions_ordered) - 1) / 2
    if 'epitope' in current_region:
        label = current_region.replace('epitope-', '')
        ax_regions.text(mid_point, -0.5, label, ha='center', va='bottom',
                       fontsize=12, fontweight='bold')

ax_regions.spines['top'].set_visible(False)
ax_regions.spines['right'].set_visible(False)
ax_regions.spines['bottom'].set_visible(False)
ax_regions.spines['left'].set_visible(False)

# Bottom panel: Heatmap
ax_heatmap = fig.add_subplot(gs[1])

# Create custom colormap (gold to white to blue)
colors_list = ['#D4A050', '#F5E6C8', '#6B9AB8']
n_bins = 100
cmap = mcolors.LinearSegmentedColormap.from_list('custom', colors_list, N=n_bins)

# Plot heatmap
vmin = -5
vmax = 2
im = ax_heatmap.imshow(heatmap_matrix.values, aspect='auto', cmap=cmap,
                       vmin=vmin, vmax=vmax, interpolation='nearest')

# Set ticks
ax_heatmap.set_xticks(range(len(heatmap_matrix.columns)))
ax_heatmap.set_xticklabels(heatmap_matrix.columns, fontsize=6, rotation=90)
ax_heatmap.set_yticks(range(len(amino_acids)))
ax_heatmap.set_yticklabels(amino_acids, fontsize=10)

# Labels
ax_heatmap.set_xlabel('Site (H7 numbering)', fontsize=14, fontweight='bold')
ax_heatmap.set_ylabel('Mutant amino acid', fontsize=14, fontweight='bold')

# Add colorbar
cbar = plt.colorbar(im, ax=ax_heatmap, fraction=0.02, pad=0.01)
cbar.set_label('Effect on cell entry', fontsize=12, fontweight='bold')

# Add title
fig.suptitle('Functional effects of mutations on HA-mediated entry into 293T 2,6 cells',
             fontsize=16, fontweight='bold', y=0.98)

# Save figure
plt.savefig('figures/full_heatmap_2-6_entry.png', dpi=300, bbox_inches='tight')
print("Heatmap saved to figures/full_heatmap_2-6_entry.png")

plt.show()
