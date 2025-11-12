"""
Create a wrapped heatmap of H7 cell entry effects on 2,6 cells with region annotations.
The heatmap will be split into 4 rows for better visualization.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors
from matplotlib.colors import to_rgb
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

# Split data into 4 rows
n_rows = 4
total_sites = len(heatmap_matrix.columns)
sites_per_row = int(np.ceil(total_sites / n_rows))

# Create the figure with 4 rows
fig = plt.figure(figsize=(20, 16))
# Each row has a small region annotation panel and a larger heatmap panel
gs = fig.add_gridspec(n_rows * 2, 1, height_ratios=[0.3, 5] * n_rows, hspace=0.05)

# Create custom colormap (gold to white to blue)
colors_list = ['#D4A050', '#F5E6C8', '#6B9AB8']
n_bins = 100
cmap = mcolors.LinearSegmentedColormap.from_list('custom', colors_list, N=n_bins)

vmin = -5
vmax = 2

# Plot each row
for row_idx in range(n_rows):
    start_idx = row_idx * sites_per_row
    end_idx = min((row_idx + 1) * sites_per_row, total_sites)

    # Get data for this row
    row_sites = heatmap_matrix.columns[start_idx:end_idx]
    row_data = heatmap_matrix[row_sites]
    row_regions = regions_ordered[start_idx:end_idx]

    # Region annotation panel
    ax_regions = fig.add_subplot(gs[row_idx * 2])
    region_bar = np.array([[to_rgb(region_colors.get(r, '#FFFFFF')) for r in row_regions]])
    ax_regions.imshow(region_bar, aspect='auto', interpolation='nearest')
    ax_regions.set_xticks([])
    ax_regions.set_yticks([])

    # Add region labels at the top
    current_region = None
    region_start = 0
    for i, region in enumerate(row_regions):
        if region != current_region:
            if current_region is not None:
                # Add label for previous region
                mid_point = (region_start + i - 1) / 2
                if 'epitope' in current_region:
                    label = current_region.replace('epitope-', '')
                    ax_regions.text(mid_point, -0.5, label, ha='center', va='bottom',
                                   fontsize=14, fontweight='bold')
            current_region = region
            region_start = i

    # Add last region label
    if current_region is not None and region_start < len(row_regions):
        mid_point = (region_start + len(row_regions) - 1) / 2
        if 'epitope' in current_region:
            label = current_region.replace('epitope-', '')
            ax_regions.text(mid_point, -0.5, label, ha='center', va='bottom',
                           fontsize=14, fontweight='bold')

    ax_regions.spines['top'].set_visible(False)
    ax_regions.spines['right'].set_visible(False)
    ax_regions.spines['bottom'].set_visible(False)
    ax_regions.spines['left'].set_visible(False)

    # Heatmap panel
    ax_heatmap = fig.add_subplot(gs[row_idx * 2 + 1])

    im = ax_heatmap.imshow(row_data.values, aspect='auto', cmap=cmap,
                           vmin=vmin, vmax=vmax, interpolation='nearest')

    # Set ticks
    ax_heatmap.set_xticks(range(len(row_sites)))
    ax_heatmap.set_xticklabels(row_sites, fontsize=8, rotation=90)
    ax_heatmap.set_yticks(range(len(amino_acids)))
    ax_heatmap.set_yticklabels(amino_acids, fontsize=12)

    # Only show x-axis label on bottom row
    if row_idx == n_rows - 1:
        ax_heatmap.set_xlabel('Site (H7 numbering)', fontsize=14, fontweight='bold')

    # Y-axis label
    ax_heatmap.set_ylabel('Mutant', fontsize=12, fontweight='bold')

# Add colorbar on the right side
cbar = plt.colorbar(im, ax=fig.axes, fraction=0.01, pad=0.01)
cbar.set_label('Effect on cell entry', fontsize=14, fontweight='bold')

# Add title
fig.suptitle('Functional effects of mutations on HA-mediated entry into 293T 2,6 cells',
             fontsize=18, fontweight='bold', y=0.995)

# Save figure
plt.savefig('figures/wrapped_heatmap_2-6_entry.png', dpi=300, bbox_inches='tight')
print("Wrapped heatmap saved to figures/wrapped_heatmap_2-6_entry.png")

plt.show()
