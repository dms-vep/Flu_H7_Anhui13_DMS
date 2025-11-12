"""
Create a wrapped heatmap of H7 cell entry effects on 2,6 cells with region annotations.
Styled to match the example heatmap with grid lines and clean layout.
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

# Merge and filter data - include sequential_site for proper numbering
func_data_for_heatmap = pd.merge(
    func_data,
    site_map[['reference_site', 'sequential_site', 'region', 'rbs_region']],
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
# Rows = amino acids, Columns = sites (using sequential_site for ordering)
amino_acids = list('ACDEFGHIKLMNPQRSTVWY')

# First, create a mapping of reference_site to sequential_site
site_to_seq = site_map.set_index('reference_site')['sequential_site'].to_dict()

# Add sequential site to the heatmap data for proper sorting
heatmap_data['sequential_site'] = heatmap_data['site'].map(site_to_seq)

# Pivot using reference site but sort by sequential site
heatmap_matrix = heatmap_data.pivot_table(
    index='mutant',
    columns='site',
    values='effect',
    aggfunc='first'
)

# Reindex to ensure all amino acids are present and in order
heatmap_matrix = heatmap_matrix.reindex(amino_acids)

# Sort columns by sequential site number
seq_site_order = heatmap_data[['site', 'sequential_site']].drop_duplicates().sort_values('sequential_site')['site'].values
heatmap_matrix = heatmap_matrix[[col for col in seq_site_order if col in heatmap_matrix.columns]]

# Create mappings for annotations
site_regions = site_map.set_index('reference_site')['region'].to_dict()
site_sequential = site_map.set_index('reference_site')['sequential_site'].to_dict()
site_wildtype = site_map.set_index('reference_site')['sequential_wt'].to_dict()
regions_ordered = [site_regions.get(site, '') for site in heatmap_matrix.columns]
sequential_sites = [site_sequential.get(site, site) for site in heatmap_matrix.columns]
wildtype_aas = [site_wildtype.get(site, '') for site in heatmap_matrix.columns]

# Define region colors
region_colors = {
    'epitope-A': '#E8B44F',
    'epitope-B': '#E8B44F',
    'epitope-C': '#E8B44F',
    'epitope-D': '#E8B44F',
    'epitope-E': '#E8B44F',
    'HA1': '#D3D3D3',
    'HA2': '#A0A0A0',
    'cleavage site': '#90C890',  # Light green for better visibility of black text
}

# Split data into 4 rows
n_rows = 4
total_sites = len(heatmap_matrix.columns)
sites_per_row = int(np.ceil(total_sites / n_rows))

# Create the figure with 4 rows
fig = plt.figure(figsize=(22, 18))
# Each row has multiple panels: region bars and heatmap
gs = fig.add_gridspec(n_rows, 1, hspace=0.4)

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
    row_sequential = sequential_sites[start_idx:end_idx]
    row_wildtype = wildtype_aas[start_idx:end_idx]

    # Create a nested gridspec for this row
    inner_gs = gs[row_idx].subgridspec(3, 1, height_ratios=[0.3, 0.3, 10], hspace=0.05)

    # Top annotation bar (epitope regions colored)
    ax_epitope = fig.add_subplot(inner_gs[0])
    epitope_bar = []
    for r in row_regions:
        if 'epitope' in r:
            epitope_bar.append(to_rgb(region_colors[r]))
        else:
            epitope_bar.append((1, 1, 1))  # White for non-epitope
    epitope_bar = np.array([epitope_bar])
    ax_epitope.imshow(epitope_bar, aspect='auto', interpolation='nearest')
    ax_epitope.set_xticks([])
    ax_epitope.set_yticks([])

    # Add epitope labels - all regions labeled with larger font
    # Collect all epitope regions first to check for overlaps
    epitope_regions = []
    current_region = None
    region_start = 0
    for i, region in enumerate(row_regions):
        if 'epitope' in region and region != current_region:
            if current_region is not None and 'epitope' in current_region:
                region_width = i - region_start
                mid_point = (region_start + i - 1) / 2
                label = current_region.replace('epitope-', '')
                epitope_regions.append((mid_point, label, region_width))
            current_region = region
            region_start = i
        elif 'epitope' not in region:
            if current_region is not None and 'epitope' in current_region:
                region_width = i - region_start
                mid_point = (region_start + i - 1) / 2
                label = current_region.replace('epitope-', '')
                epitope_regions.append((mid_point, label, region_width))
            current_region = None

    # Add last epitope label if needed
    if current_region is not None and 'epitope' in current_region:
        region_width = len(row_regions) - region_start
        mid_point = (region_start + len(row_regions) - 1) / 2
        label = current_region.replace('epitope-', '')
        epitope_regions.append((mid_point, label, region_width))

    # Draw labels with larger font, checking for overlaps
    fontsize = 14  # Larger base font size
    for i, (mid_point, label, region_width) in enumerate(epitope_regions):
        # Check distance to previous label to avoid overlap
        if i > 0:
            prev_mid = epitope_regions[i-1][0]
            distance = mid_point - prev_mid
            if distance < 20:  # If too close, reduce font size slightly
                fontsize = 12
            else:
                fontsize = 14
        ax_epitope.text(mid_point, 0.2, label, ha='center', va='center',
                       fontsize=fontsize, fontweight='bold')

    for spine in ax_epitope.spines.values():
        spine.set_visible(False)

    # Second annotation bar (HA1/HA2 regions)
    ax_ha = fig.add_subplot(inner_gs[1])
    ha_bar = []
    for r in row_regions:
        if 'HA' in r or 'cleavage' in r:
            ha_bar.append(to_rgb(region_colors.get(r, '#FFFFFF')))
        else:
            ha_bar.append((1, 1, 1))
    ha_bar = np.array([ha_bar])
    ax_ha.imshow(ha_bar, aspect='auto', interpolation='nearest')
    ax_ha.set_xticks([])
    ax_ha.set_yticks([])

    # Add HA1/HA2 labels - collect regions first to check overlaps
    ha_regions = []
    current_region = None
    region_start = 0
    for i, region in enumerate(row_regions):
        if region in ['HA1', 'HA2', 'cleavage site'] and region != current_region:
            if current_region is not None:
                region_width = i - region_start
                mid_point = (region_start + i - 1) / 2
                label = current_region
                if current_region == 'cleavage site':
                    label = 'PBCS'
                ha_regions.append((mid_point, label, region_width))
            current_region = region
            region_start = i
        elif region not in ['HA1', 'HA2', 'cleavage site']:
            if current_region is not None:
                region_width = i - region_start
                mid_point = (region_start + i - 1) / 2
                label = current_region
                if current_region == 'cleavage site':
                    label = 'PBCS'
                ha_regions.append((mid_point, label, region_width))
            current_region = None

    if current_region is not None:
        region_width = len(row_regions) - region_start
        mid_point = (region_start + len(row_regions) - 1) / 2
        label = current_region
        if current_region == 'cleavage site':
            label = 'PBCS'
        ha_regions.append((mid_point, label, region_width))

    # Draw labels with larger font, checking for overlaps
    fontsize_base = 14  # Larger base font size
    for i, (mid_point, label, region_width) in enumerate(ha_regions):
        # Check distance to previous label to avoid overlap
        if i > 0:
            prev_mid = ha_regions[i-1][0]
            distance = mid_point - prev_mid
            if distance < 25:  # If too close, reduce font size
                fontsize = 11
            else:
                fontsize = fontsize_base
        else:
            fontsize = fontsize_base
        ax_ha.text(mid_point, 0.2, label, ha='center', va='center',
                  fontsize=fontsize, fontweight='normal', color='black')

    for spine in ax_ha.spines.values():
        spine.set_visible(False)

    # Heatmap panel
    ax_heatmap = fig.add_subplot(inner_gs[2])

    im = ax_heatmap.imshow(row_data.values, aspect='auto', cmap=cmap,
                           vmin=vmin, vmax=vmax, interpolation='nearest')

    # Add grid lines
    ax_heatmap.set_xticks(np.arange(len(row_sites)) - 0.5, minor=True)
    ax_heatmap.set_yticks(np.arange(len(amino_acids)) - 0.5, minor=True)
    ax_heatmap.grid(which='minor', color='gray', linestyle='-', linewidth=0.5)
    ax_heatmap.tick_params(which='minor', size=0)

    # Set ticks - only label every 20th site
    tick_positions = []
    tick_labels = []
    for i, seq_site in enumerate(row_sequential):
        if seq_site % 20 == 0:  # Label every 20th site
            tick_positions.append(i)
            tick_labels.append(str(seq_site))

    ax_heatmap.set_xticks(tick_positions)
    ax_heatmap.set_xticklabels(tick_labels, fontsize=9)
    ax_heatmap.set_yticks(range(len(amino_acids)))
    ax_heatmap.set_yticklabels(amino_acids, fontsize=10, family='monospace')

    # Mark wildtype positions with X
    for col_idx, wt_aa in enumerate(row_wildtype):
        if wt_aa in amino_acids:
            row_idx_wt = amino_acids.index(wt_aa)
            ax_heatmap.text(col_idx, row_idx_wt, 'X', ha='center', va='center',
                          fontsize=6, fontweight='bold', color='black')

    # Position x-axis labels at bottom
    ax_heatmap.xaxis.tick_bottom()

    # Y-axis label
    if row_idx == 0:
        ax_heatmap.set_ylabel('Mutant amino acid', fontsize=12, fontweight='bold')
    else:
        ax_heatmap.set_ylabel('')

    # Remove top and right spines
    ax_heatmap.spines['top'].set_visible(False)
    ax_heatmap.spines['right'].set_visible(False)
    ax_heatmap.spines['bottom'].set_visible(True)
    ax_heatmap.spines['left'].set_visible(True)

# Add colorbar on the right side
cbar_ax = fig.add_axes([0.92, 0.15, 0.015, 0.3])
cbar = fig.colorbar(im, cax=cbar_ax)
cbar.set_label('Effect on cell entry', fontsize=12, fontweight='bold')

# Add title
fig.suptitle('Functional effects of mutations on HA-mediated entry into 293T 2,6 cells',
             fontsize=16, fontweight='bold', y=0.98)

# Save figure
plt.savefig('figures/styled_heatmap_2-6_entry.png', dpi=300, bbox_inches='tight', facecolor='white')
print("Styled heatmap saved to figures/styled_heatmap_2-6_entry.png")

plt.show()
