#!/usr/bin/env python3
"""
Verify that the regenerated 293_2-6_entry_func_effects.csv matches
the existing defattr files and numbering table.
"""

import pandas as pd
import numpy as np

print("=" * 80)
print("VERIFICATION: Regenerated func_effects data vs existing files")
print("=" * 80)

# 1. Read regenerated func_effects (reference numbering)
print("\n1. Reading regenerated func_effects file...")
func_regen = pd.read_csv('results/func_effects/averages/293_2-6_entry_func_effects.csv')
print(f"   Rows: {len(func_regen)}")
print(f"   Unique sites: {func_regen['site'].nunique()}")

# Calculate averages from regenerated file
site_avg_regen = func_regen.groupby('site')['effect'].mean().reset_index()
site_avg_regen.columns = ['reference_site', 'avg_effect_regen']
site_avg_regen['reference_site'] = site_avg_regen['reference_site'].astype(str)

# 2. Read existing verification table
print("\n2. Reading existing verification table...")
table = pd.read_csv('site_numbering_verification_table.csv')
table['reference_site'] = table['reference_site'].astype(str)
print(f"   Rows: {len(table)}")

# 3. Merge and compare
print("\n3. Comparing values...")
comparison = pd.merge(
    table[['sequential_site', 'reference_site', 'sequential_wt', 'region', 'avg_func_effect']],
    site_avg_regen,
    on='reference_site',
    how='left'
)

# Calculate difference
comparison['difference'] = comparison['avg_func_effect'] - comparison['avg_effect_regen']

# Check for discrepancies
max_diff = comparison['difference'].abs().max()
print(f"   Maximum difference: {max_diff:.10f}")

if max_diff < 1e-6:
    print("   ✅ ALL VALUES MATCH (within rounding error)")
else:
    print(f"   ⚠️  Some differences found > 1e-6")
    discrepancies = comparison[comparison['difference'].abs() > 1e-6]
    print(f"   Sites with discrepancies: {len(discrepancies)}")
    print("\n   First 10 discrepancies:")
    print(discrepancies[['sequential_site', 'reference_site', 'avg_func_effect',
                         'avg_effect_regen', 'difference']].head(10))

# 4. Check previously missing sites
print("\n4. Checking previously missing sites (147, 253, 331, 332)...")
missing_check = comparison[comparison['sequential_site'].isin([147, 253, 331, 332])]
print(missing_check[['sequential_site', 'reference_site', 'sequential_wt',
                     'avg_func_effect', 'avg_effect_regen']].to_string(index=False))

if missing_check['avg_effect_regen'].isna().any():
    print("   ⚠️  Some previously missing sites are STILL missing!")
else:
    print("   ✅ All previously missing sites are now present")

# 5. Check specific test cases
print("\n5. Testing specific sites:")
test_cases = [
    (1, '11', 'First site'),
    (147, '158', 'Previously missing'),
    (253, '261a', 'Previously missing'),
    (317, '326', 'Last HA1'),
    (323, '331', 'First HA2'),
    (331, '339', 'Previously missing'),
    (332, '340', 'Previously missing'),
    (506, '514', 'Last site')
]

for seq_site, ref_site, description in test_cases:
    row = comparison[comparison['sequential_site'] == seq_site]
    if not row.empty:
        old_val = row['avg_func_effect'].values[0]
        new_val = row['avg_effect_regen'].values[0]
        diff = abs(old_val - new_val) if pd.notna(new_val) else float('nan')
        status = "✅" if pd.notna(new_val) and diff < 1e-6 else "⚠️"
        print(f"   {status} Seq {seq_site} (Ref {ref_site}) - {description}")
        print(f"      Old: {old_val:.4f}, New: {new_val:.4f}, Diff: {diff:.10f}")

# 6. Summary
print("\n" + "=" * 80)
print("SUMMARY")
print("=" * 80)
print(f"Total sites compared: {len(comparison)}")
print(f"Sites matching: {(comparison['difference'].abs() < 1e-6).sum()}")
print(f"Sites with discrepancies: {(comparison['difference'].abs() >= 1e-6).sum()}")
print(f"Missing in regenerated file: {comparison['avg_effect_regen'].isna().sum()}")

if max_diff < 1e-6 and not comparison['avg_effect_regen'].isna().any():
    print("\n✅ VERIFICATION PASSED: All data matches perfectly!")
else:
    print("\n⚠️  VERIFICATION FAILED: Some issues detected")
