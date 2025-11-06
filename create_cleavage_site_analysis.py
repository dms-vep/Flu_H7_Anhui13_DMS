"""
Detailed analysis of cleavage site mutations showing positive vs negative effects.
Creates a comprehensive plot highlighting the unusual positive effects at K328.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Read the data
func_data = pd.read_csv('results/func_effects/averages/293_mix_entry_func_effects.csv')
site_map = pd.read_csv('data/site_numbering_map.csv')

# Merge and filter
func_data_filtered = pd.merge(
    func_data,
    site_map[['reference_site', 'sequential_site', 'region']],
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

# Get cleavage site mutations
cleavage_muts = func_data_filtered[func_data_filtered['region'] == 'cleavage site'].copy()

# Create figure with multiple panels
fig = plt.figure(figsize=(18, 10))
gs = fig.add_gridspec(2, 3, hspace=0.3, wspace=0.3)

# Panel 1: Distribution of effects by site
ax1 = fig.add_subplot(gs[0, :])
sites_order = ['327', '328', '328a', '329']
site_data = []
for site in sites_order:
    site_muts = cleavage_muts[cleavage_muts['site'] == site]
    for _, row in site_muts.iterrows():
        site_data.append({
            'site': site,
            'effect': row['effect'],
            'mutation': f"{row['wildtype']}{site}{row['mutant']}",
            'is_positive': row['effect'] > 0
        })

site_df = pd.DataFrame(site_data)

# Create violin plot with swarm overlay
parts = ax1.violinplot(
    [site_df[site_df['site'] == s]['effect'].values for s in sites_order],
    positions=range(len(sites_order)),
    showmeans=True,
    showmedians=True,
    widths=0.7
)

# Color violins
for pc in parts['bodies']:
    pc.set_facecolor('#1A391A')
    pc.set_alpha(0.3)

# Add scatter points colored by positive/negative
for i, site in enumerate(sites_order):
    site_vals = site_df[site_df['site'] == site]
    positive = site_vals[site_vals['is_positive']]
    negative = site_vals[~site_vals['is_positive']]

    # Jitter x positions
    jitter = 0.1
    pos_x = np.random.normal(i, jitter, size=len(positive))
    neg_x = np.random.normal(i, jitter, size=len(negative))

    ax1.scatter(pos_x, positive['effect'], c='#E74C3C', s=50, alpha=0.7,
                edgecolors='black', linewidth=0.5, label='Positive' if i == 0 else '')
    ax1.scatter(neg_x, negative['effect'], c='#3498DB', s=50, alpha=0.7,
                edgecolors='black', linewidth=0.5, label='Negative' if i == 0 else '')

ax1.axhline(0, color='black', linestyle='--', linewidth=1, alpha=0.5)
ax1.set_xticks(range(len(sites_order)))
ax1.set_xticklabels([f'Site {s}' for s in sites_order], fontsize=12, fontweight='bold')
ax1.set_ylabel('Effect on cell entry', fontsize=13, fontweight='bold')
ax1.set_title('Distribution of mutation effects across cleavage site positions',
              fontsize=14, fontweight='bold')
ax1.legend(loc='upper right', fontsize=11)
ax1.grid(True, alpha=0.3, axis='y')

# Panel 2: Site 328 detailed view (most positive effects)
ax2 = fig.add_subplot(gs[1, 0])
k328_muts = cleavage_muts[cleavage_muts['site'] == '328'].sort_values('effect', ascending=True)
colors = ['#E74C3C' if x > 0 else '#3498DB' for x in k328_muts['effect']]
y_pos = np.arange(len(k328_muts))

ax2.barh(y_pos, k328_muts['effect'], color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
ax2.set_yticks(y_pos)
ax2.set_yticklabels([f"{row['wildtype']}{row['site']}{row['mutant']}"
                      for _, row in k328_muts.iterrows()], fontsize=9)
ax2.axvline(0, color='black', linestyle='-', linewidth=1)
ax2.set_xlabel('Effect on cell entry', fontsize=11, fontweight='bold')
ax2.set_title('Site 328 (K) mutations\nMost positive effects', fontsize=12, fontweight='bold')
ax2.grid(True, alpha=0.3, axis='x')

# Panel 3: Site 329 detailed view (mostly negative)
ax3 = fig.add_subplot(gs[1, 1])
r329_muts = cleavage_muts[cleavage_muts['site'] == '329'].sort_values('effect', ascending=True).head(15)
colors = ['#E74C3C' if x > 0 else '#3498DB' for x in r329_muts['effect']]
y_pos = np.arange(len(r329_muts))

ax3.barh(y_pos, r329_muts['effect'], color=colors, alpha=0.7, edgecolor='black', linewidth=0.5)
ax3.set_yticks(y_pos)
ax3.set_yticklabels([f"{row['wildtype']}{row['site']}{row['mutant']}"
                      for _, row in r329_muts.iterrows()], fontsize=9)
ax3.axvline(0, color='black', linestyle='-', linewidth=1)
ax3.set_xlabel('Effect on cell entry', fontsize=11, fontweight='bold')
ax3.set_title('Site 329 (R) mutations\n(Top 15 most deleterious)', fontsize=12, fontweight='bold')
ax3.grid(True, alpha=0.3, axis='x')

# Panel 4: Summary statistics
ax4 = fig.add_subplot(gs[1, 2])
ax4.axis('off')

# Calculate statistics per site
stats_text = "Summary Statistics:\n\n"
for site in sites_order:
    site_muts = cleavage_muts[cleavage_muts['site'] == site]
    wt = site_muts.iloc[0]['wildtype'] if len(site_muts) > 0 else '?'
    n_total = len(site_muts)
    n_pos = len(site_muts[site_muts['effect'] > 0])
    n_neg = len(site_muts[site_muts['effect'] < 0])
    mean_eff = site_muts['effect'].mean()
    median_eff = site_muts['effect'].median()

    stats_text += f"Site {site} ({wt}):\n"
    stats_text += f"  Total: {n_total} mutations\n"
    stats_text += f"  Positive: {n_pos} ({100*n_pos/n_total:.1f}%)\n"
    stats_text += f"  Negative: {n_neg} ({100*n_neg/n_total:.1f}%)\n"
    stats_text += f"  Mean: {mean_eff:.3f}\n"
    stats_text += f"  Median: {median_eff:.3f}\n\n"

ax4.text(0.1, 0.9, stats_text, transform=ax4.transAxes,
         fontsize=10, verticalalignment='top', fontfamily='monospace',
         bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.3))

# Overall title
fig.suptitle('Cleavage Site Mutation Analysis: 293T Mix Cell Entry\n' +
             'Polybasic Cleavage Site (PBCS): P327-K328-G328a-R329',
             fontsize=16, fontweight='bold', y=0.98)

plt.savefig('figures/cleavage_site_detailed_analysis.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/cleavage_site_detailed_analysis.svg', bbox_inches='tight')
print("Saved detailed cleavage site analysis plots")
print(f"figures/cleavage_site_detailed_analysis.png")
print(f"figures/cleavage_site_detailed_analysis.svg")
print()

# Create summary table
print("\n" + "="*80)
print("CLEAVAGE SITE MUTATION SUMMARY")
print("="*80)
print(f"\nTotal cleavage site mutations: {len(cleavage_muts)}")
print(f"Mutations with positive effects: {len(cleavage_muts[cleavage_muts['effect'] > 0])} ({100*len(cleavage_muts[cleavage_muts['effect'] > 0])/len(cleavage_muts):.1f}%)")
print(f"Mutations with negative effects: {len(cleavage_muts[cleavage_muts['effect'] < 0])} ({100*len(cleavage_muts[cleavage_muts['effect'] < 0])/len(cleavage_muts):.1f}%)")
print(f"\nMean effect: {cleavage_muts['effect'].mean():.3f}")
print(f"Median effect: {cleavage_muts['effect'].median():.3f}")

print("\n" + "-"*80)
print("TOP 10 POSITIVE EFFECTS (Enhanced entry):")
print("-"*80)
for idx, row in cleavage_muts[cleavage_muts['effect'] > 0].sort_values('effect', ascending=False).head(10).iterrows():
    print(f"{row['wildtype']}{row['site']}{row['mutant']:>3s}: {row['effect']:7.3f}  (std={row['effect_std']:.3f}, n={row['times_seen']:.1f})")

print("\n" + "-"*80)
print("TOP 10 NEGATIVE EFFECTS (Reduced entry):")
print("-"*80)
for idx, row in cleavage_muts.sort_values('effect').head(10).iterrows():
    print(f"{row['wildtype']}{row['site']}{row['mutant']:>3s}: {row['effect']:7.3f}  (std={row['effect_std']:.3f}, n={row['times_seen']:.1f})")

print("\n" + "="*80)
print("KEY BIOLOGICAL FINDING:")
print("="*80)
print("""
Site 328 (K - Lysine):
  - Removing the positive charge (K→D, K→E, K→A, K→G, etc.) ENHANCES entry
  - This suggests that reducing positive charge at K328 may improve cleavage
    efficiency by furin-like proteases

Site 329 (R - Arginine):
  - Almost ALL mutations are highly deleterious
  - R329 is CRITICAL for furin recognition and cleavage
  - This is the key determinant for polybasic cleavage site function

Interpretation:
  The R329 is absolutely required, but K328 may be sub-optimal.
  Natural selection may favor K328 for reasons other than cleavage efficiency.
""")
print("="*80)
