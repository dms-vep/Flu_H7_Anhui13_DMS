"""
Create summary plot of cell entry effects across all three cell types.
Apply consistent filters (times_seen >= 2, effect_std <= 2) and save filtered data.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Define the three datasets
datasets = {
    '293T 2-3': 'results/func_effects/averages/293_2-3_entry_func_effects.csv',
    '293T 2-6': 'results/func_effects/averages/293_2-6_entry_func_effects.csv',
    '293T mix': 'results/func_effects/averages/293_mix_entry_func_effects.csv',
}

# Read and filter all datasets
filtered_data = {}
all_mutations = set()

print("Applying filters: times_seen >= 2 AND effect_std <= 2\n")

for name, path in datasets.items():
    # Read data
    df = pd.read_csv(path)

    # Filter
    df_filtered = df[
        (df['times_seen'] >= 2) &
        (df['effect_std'] <= 2) &
        (df['mutant'] != '*') &  # Remove stop codons
        (df['mutant'] != df['wildtype'])  # Remove wildtype
    ].copy()

    # Create mutation identifier
    df_filtered['mutation'] = df_filtered['wildtype'] + df_filtered['site'].astype(str) + df_filtered['mutant']

    # Store filtered data
    filtered_data[name] = df_filtered
    all_mutations.update(df_filtered['mutation'].tolist())

    # Print stats
    print(f"{name}:")
    print(f"  Total mutations before filtering: {len(df)}")
    print(f"  Mutations after filtering: {len(df_filtered)}")
    print(f"  Filtered out: {len(df) - len(df_filtered)} ({100*(len(df) - len(df_filtered))/len(df):.1f}%)")
    print()

# Save filtered datasets
print("Saving filtered datasets...")
for name, df in filtered_data.items():
    output_name = name.replace(' ', '_').replace('T', 'T_')
    output_path = f'results/func_effects/filtered/{output_name}_filtered.csv'
    df.to_csv(output_path, index=False)
    print(f"  Saved: {output_path}")

print()

# Create combined dataset for plotting
combined_data = []
for name, df in filtered_data.items():
    temp = df[['mutation', 'site', 'wildtype', 'mutant', 'effect', 'effect_std']].copy()
    temp['cell_type'] = name
    combined_data.append(temp)

combined_df = pd.concat(combined_data, ignore_index=True)

# Save combined filtered data
combined_output = 'results/func_effects/filtered/all_cells_entry_filtered.csv'
combined_df.to_csv(combined_output, index=False)
print(f"Saved combined filtered data: {combined_output}\n")

# Create pairwise correlation matrix
fig, axes = plt.subplots(2, 2, figsize=(14, 14))

# Merge datasets for pairwise comparisons
datasets_list = list(filtered_data.items())
comparisons = [
    (0, 1, axes[0, 0]),  # 2-3 vs 2-6
    (0, 2, axes[0, 1]),  # 2-3 vs mix
    (1, 2, axes[1, 0]),  # 2-6 vs mix
]

for idx1, idx2, ax in comparisons:
    name1, df1 = datasets_list[idx1]
    name2, df2 = datasets_list[idx2]

    # Merge on mutation
    merged = pd.merge(
        df1[['mutation', 'effect']],
        df2[['mutation', 'effect']],
        on='mutation',
        suffixes=(f'_{name1}', f'_{name2}')
    )

    # Calculate correlation
    r = np.corrcoef(merged[f'effect_{name1}'], merged[f'effect_{name2}'])[0, 1]

    # Create scatter plot
    ax.scatter(merged[f'effect_{name1}'], merged[f'effect_{name2}'],
               alpha=0.3, s=10, color='steelblue')

    # Add diagonal line
    lims = [
        min(merged[f'effect_{name1}'].min(), merged[f'effect_{name2}'].min()),
        max(merged[f'effect_{name1}'].max(), merged[f'effect_{name2}'].max())
    ]
    ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=1)

    # Labels and title
    ax.set_xlabel(f'{name1} entry effect', fontsize=12, fontweight='bold')
    ax.set_ylabel(f'{name2} entry effect', fontsize=12, fontweight='bold')
    ax.set_title(f'n={len(merged):,} mutations\nPearson r={r:.3f} (p<0.001)',
                 fontsize=11, fontweight='bold')
    ax.grid(True, alpha=0.3)
    ax.axhline(0, color='gray', linewidth=0.8, alpha=0.5)
    ax.axvline(0, color='gray', linewidth=0.8, alpha=0.5)

# Remove unused subplot
axes[1, 1].remove()

# Add overall title
fig.suptitle('Cell entry functional effects correlations across cell types\n(Filtered: times_seen ≥ 2, effect_std ≤ 2)',
             fontsize=16, fontweight='bold', y=0.995)

plt.tight_layout()
plt.savefig('figures/entry_summary_correlations.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/entry_summary_correlations.svg', bbox_inches='tight')
print("Saved correlation plots: figures/entry_summary_correlations.png/.svg")

# Create distribution comparison plot
fig, axes = plt.subplots(1, 3, figsize=(18, 5))

for idx, (name, df) in enumerate(filtered_data.items()):
    ax = axes[idx]

    # Plot distribution
    ax.hist(df['effect'], bins=50, alpha=0.7, color='steelblue', edgecolor='black')
    ax.axvline(0, color='red', linestyle='--', linewidth=2, label='Neutral (0)')

    # Add stats
    mean_effect = df['effect'].mean()
    median_effect = df['effect'].median()
    ax.axvline(mean_effect, color='orange', linestyle='-', linewidth=2, label=f'Mean ({mean_effect:.2f})')

    ax.set_xlabel('Entry effect', fontsize=12, fontweight='bold')
    ax.set_ylabel('Count', fontsize=12, fontweight='bold')
    ax.set_title(f'{name}\n{len(df):,} mutations', fontsize=13, fontweight='bold')
    ax.legend(fontsize=10)
    ax.grid(True, alpha=0.3, axis='y')

fig.suptitle('Distribution of cell entry effects across cell types\n(Filtered: times_seen ≥ 2, effect_std ≤ 2)',
             fontsize=16, fontweight='bold')

plt.tight_layout()
plt.savefig('figures/entry_summary_distributions.png', dpi=300, bbox_inches='tight')
plt.savefig('figures/entry_summary_distributions.svg', bbox_inches='tight')
print("Saved distribution plots: figures/entry_summary_distributions.png/.svg")

print("\nSummary complete!")
