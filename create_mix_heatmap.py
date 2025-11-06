"""
Create a wrapped heatmap of H7 cell entry effects on 293 mix cells with region annotations.
Styled to match the example heatmap with grid lines and clean layout.
Color scale: centered at 0 (white), missing data in grey.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib import colors as mcolors
from matplotlib.colors import to_rgb
import seaborn as sns

# Read in the data
func_data = pd.read_csv('results/func_effects/averages/293_mix_entry_func_effects.csv')
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
# Order amino acids biochemically: Small, Hydrophobic, Aromatic, Polar, Basic, Acidic
amino_acids = list('GAPVILMFWYSTEDQNKHRCDE')

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

# Define region colors - use distinct colors not in heatmap (blues, greens, purples)
region_colors = {
    'epitope-A': '#9DC3E6',  # Light blue
    'epitope-B': '#9DC3E6',  # Light blue
    'epitope-C': '#9DC3E6',  # Light blue
    'epitope-D': '#9DC3E6',  # Light blue
    'epitope-E': '#9DC3E6',  # Light blue
    'HA1': '#C5E0B4',  # Light green
    'HA2': '#D9D2E9',  # Light purple
    'cleavage site': '#FFE699',  # Light yellow
}

# Split data into 3 rows
n_rows = 3
total_sites = len(heatmap_matrix.columns)
sites_per_row = int(np.ceil(total_sites / n_rows))

# Create the figure with 3 rows - wider to accommodate more sites per row, taller for more spacing
fig = plt.figure(figsize=(32, 20))
# Each row has multiple panels: region bars and heatmap
# Minimize left whitespace, reduce bottom margin for less space above colorbar
gs = fig.add_gridspec(n_rows, 1, hspace=0.25, left=0.04, right=0.96, top=0.95, bottom=0.06)

# Create custom colormap - centered at 0 (white), with grey for missing data
colors_list = ['#D4A050', '#F5E6C8', 'white', '#A4C8D8', '#6B9AB8']
n_bins = 100
cmap = mcolors.LinearSegmentedColormap.from_list('custom', colors_list, N=n_bins)
cmap.set_bad(color='#CCCCCC')  # Grey for missing/NaN values

# Determine vmin and vmax centered at 0
vmax = max(abs(heatmap_matrix.min().min()), abs(heatmap_matrix.max().max()))
vmin = -vmax

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

    # Create a nested gridspec for this row - only 3 bars now
    inner_gs = gs[row_idx].subgridspec(3, 1, height_ratios=[0.3, 0.3, 10], hspace=0.05)

    # Top annotation bar (HA1/HA2 continuous regions)
    ax_ha_continuous = fig.add_subplot(inner_gs[0])
    ha_continuous_bar = []
    for r in row_regions:
        if 'HA1' in r or 'epitope' in r:  # Group all HA1 including epitopes
            ha_continuous_bar.append(to_rgb(region_colors['HA1']))
        elif 'HA2' in r:
            ha_continuous_bar.append(to_rgb(region_colors['HA2']))
        elif 'cleavage' in r:
            ha_continuous_bar.append(to_rgb(region_colors['cleavage site']))
        else:
            ha_continuous_bar.append((1, 1, 1))
    ha_continuous_bar = np.array([ha_continuous_bar])
    ax_ha_continuous.imshow(ha_continuous_bar, aspect='auto', interpolation='nearest')
    ax_ha_continuous.set_xticks([])
    ax_ha_continuous.set_yticks([])
    ax_ha_continuous.set_frame_on(False)  # Remove frame completely

    # Add HA1/HA2/PBCS labels for continuous bar
    ha_cont_regions = []
    current_region = None
    region_start = 0
    for i, region in enumerate(row_regions):
        # Classify into HA1, HA2, or PBCS
        if ('HA1' in region or 'epitope' in region) and current_region != 'HA1':
            if current_region is not None:
                mid_point = (region_start + i - 1) / 2
                ha_cont_regions.append((mid_point, current_region))
            current_region = 'HA1'
            region_start = i
        elif 'HA2' in region and current_region != 'HA2':
            if current_region is not None:
                mid_point = (region_start + i - 1) / 2
                ha_cont_regions.append((mid_point, current_region))
            current_region = 'HA2'
            region_start = i
        elif 'cleavage' in region and current_region != 'PBCS':
            if current_region is not None:
                mid_point = (region_start + i - 1) / 2
                ha_cont_regions.append((mid_point, current_region))
            current_region = 'PBCS'
            region_start = i

    if current_region is not None:
        mid_point = (region_start + len(row_regions) - 1) / 2
        ha_cont_regions.append((mid_point, current_region))

    # Draw continuous HA labels
    for mid_point, label in ha_cont_regions:
        ax_ha_continuous.text(mid_point, 0.2, label, ha='center', va='center',
                             fontsize=19, fontweight='normal', color='black')

    for spine in ax_ha_continuous.spines.values():
        spine.set_visible(False)
    # Also turn off axis lines
    ax_ha_continuous.axhline(visible=False)
    ax_ha_continuous.axvline(visible=False)

    # Second annotation bar (epitope regions colored)
    ax_epitope = fig.add_subplot(inner_gs[1])
    epitope_bar = []
    for r in row_regions:
        if 'epitope' in r:
            # Use the same light blue color for all epitopes
            epitope_bar.append(to_rgb('#9DC3E6'))
        else:
            epitope_bar.append((1, 1, 1))  # White for non-epitope
    epitope_bar = np.array([epitope_bar])
    ax_epitope.imshow(epitope_bar, aspect='auto', interpolation='nearest')
    ax_epitope.set_xticks([])
    ax_epitope.set_yticks([])
    ax_epitope.set_frame_on(False)  # Remove frame completely

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
    fontsize = 19  # Larger base font size (16 * 1.2)
    for i, (mid_point, label, region_width) in enumerate(epitope_regions):
        # Check distance to previous label to avoid overlap
        if i > 0:
            prev_mid = epitope_regions[i-1][0]
            distance = mid_point - prev_mid
            if distance < 20:  # If too close, reduce font size slightly
                fontsize = 17
            else:
                fontsize = 19
        ax_epitope.text(mid_point, 0.2, label, ha='center', va='center',
                       fontsize=fontsize, fontweight='bold')

    for spine in ax_epitope.spines.values():
        spine.set_visible(False)
    # Also turn off axis lines
    ax_epitope.axhline(visible=False)
    ax_epitope.axvline(visible=False)

    # Heatmap panel (now third row, not fourth)
    ax_heatmap = fig.add_subplot(inner_gs[2])
    ax_heatmap.set_frame_on(False)  # Remove frame completely from heatmap too

    im = ax_heatmap.imshow(row_data.values, aspect='auto', cmap=cmap,
                           vmin=vmin, vmax=vmax, interpolation='nearest')

    # Add grid lines
    ax_heatmap.set_xticks(np.arange(len(row_sites)) - 0.5, minor=True)
    ax_heatmap.set_yticks(np.arange(len(amino_acids)) - 0.5, minor=True)
    ax_heatmap.grid(which='minor', color='gray', linestyle='-', linewidth=0.5)
    ax_heatmap.tick_params(which='minor', size=0)
    # Turn off any axis lines
    ax_heatmap.axhline(y=-0.5, visible=False)
    ax_heatmap.axhline(y=len(amino_acids)-0.5, visible=False)

    # Set ticks - only label every 20th site with larger font
    tick_positions = []
    tick_labels = []
    for i, seq_site in enumerate(row_sequential):
        if seq_site % 20 == 0:  # Label every 20th site
            tick_positions.append(i)
            tick_labels.append(str(seq_site))

    # Show x-axis ticks and labels on all rows
    ax_heatmap.set_xticks(tick_positions)
    ax_heatmap.set_xticklabels(tick_labels, fontsize=19)  # Larger font for numbers (16 * 1.2)
    ax_heatmap.xaxis.tick_bottom()

    ax_heatmap.set_yticks(range(len(amino_acids)))
    ax_heatmap.set_yticklabels(amino_acids, fontsize=19, family='monospace')  # Larger font (16 * 1.2)

    # Mark wildtype positions with X
    for col_idx, wt_aa in enumerate(row_wildtype):
        if wt_aa in amino_acids:
            row_idx_wt = amino_acids.index(wt_aa)
            ax_heatmap.text(col_idx, row_idx_wt, 'X', ha='center', va='center',
                          fontsize=8, fontweight='bold', color='black')

    # Y-axis label - centered across all rows
    ax_heatmap.set_ylabel('')

    # Remove all spines from heatmap to avoid dark bars
    for spine in ax_heatmap.spines.values():
        spine.set_visible(False)
        spine.set_linewidth(0)

    # Completely turn off axis lines and ticks
    ax_heatmap.tick_params(axis='both', which='both', length=0, width=0,
                           bottom=False, top=False, left=False, right=False,
                           labelbottom=True, labelleft=True)
    ax_heatmap.xaxis.set_ticks_position('none')
    ax_heatmap.yaxis.set_ticks_position('none')

    # Disable any drawn axis lines
    for side in ['top', 'bottom', 'left', 'right']:
        ax_heatmap.spines[side].set_color('none')

# Add colorbar at the bottom - aligned with right side of heatmap
cbar_ax = fig.add_axes([0.76, 0.015, 0.2, 0.015])  # [left, bottom, width, height]
cbar = fig.colorbar(im, cax=cbar_ax, orientation='horizontal')
cbar.set_label('Effect on cell entry', fontsize=17, fontweight='bold', labelpad=8)
cbar.ax.tick_params(labelsize=14)

# Add centered Y-axis label for entire figure - positioned at 0.02
fig.text(0.02, 0.5, 'Mutant amino acid', va='center', rotation='vertical',
         fontsize=22, fontweight='bold')

# Add title with more whitespace at top
fig.suptitle('Functional effects of mutations on HA-mediated entry into 293T mix cells',
             fontsize=24, fontweight='bold', y=0.98)

# Save figure with minimal whitespace
plt.savefig('figures/mix_heatmap_entry.png', dpi=300, bbox_inches='tight',
            facecolor='white', pad_inches=0.1)
print("Mix heatmap saved to figures/mix_heatmap_entry.png")

# Save as SVG for editing in Illustrator
plt.savefig('figures/mix_heatmap_entry.svg', bbox_inches='tight',
            facecolor='white', pad_inches=0.1)
print("Mix heatmap saved to figures/mix_heatmap_entry.svg")

plt.show()
