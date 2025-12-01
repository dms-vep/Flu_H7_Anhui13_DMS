ChimeraX Attribute Files for H7 Anhui13 Cell Entry Effects
===========================================================

Files:
------
1. chimera_293_mix_entry.defattr
   - Attribute: cell_entry_effect (numeric values)
   - Contains mean cell entry effect for each residue
   - Range: -4.789 to 0.699
   - Applied to all 3 protomers in the trimer
   
2. chimera_293_mix_entry_color.defattr
   - Attribute: cell_entry_color (RGB colors)
   - Red-blue color scheme matching wrapped heatmap
   - Red = negative effects (reduced entry)
   - White = neutral effects
   - Blue = positive effects (enhanced entry)
   - Fixed range: -5.014 to 1.426
   - Applied to all 3 protomers in the trimer

Data Source:
-----------
- Source: 293_mix_entry_func_effects.csv
- Processing: Averaged effect values across all mutations per site
- Total sites: 506 residues
- Numbering: Sequential PDB numbering matching 6ii9_JA_renumbered.pdb

Structure:
----------
6ii9_JA_renumbered.pdb contains an HA trimer with 6 chains:
- Chains A, C, E: HA1 (317 residues each) - 3 copies
- Chains B, D, F: HA2 (170 residues each) - 3 copies

Attributes are applied to ALL protomers so the entire trimer is colored.

Usage in ChimeraX:
------------------
1. Open structure:
   open analysis/chimera/6ii9_JA_renumbered.pdb

2. Load attribute files:
   open chimera_293_mix_entry.defattr
   open chimera_293_mix_entry_color.defattr

3. Color entire trimer by attribute:
   color byattribute cell_entry_color

4. Color just one protomer (chains A and B):
   color /A,B byattribute cell_entry_color
   color /C-F tan

5. Display values (optional):
   label /A,B text {0.cell_entry_effect:.2f}

6. Use in selection/analysis:
   select ::cell_entry_effect<-2
   select ::cell_entry_effect>0
   select /A:90-150::cell_entry_effect<-3

Format:
-------
- Standard ChimeraX .defattr format
- Residue-level attributes
- Atom specification: /chain:residue (e.g., /A:90, /B:322)
- Match mode: any (allows multiple matches per residue number)
- Tab-delimited assignment lines
- Total assignments: 1518 (506 unique sites × 3 protomers)

Color Scheme:
-------------
Matches the wrapped heatmap red-blue color scheme:
- Values ≤ -5.014: Bright red
- Values = -2.5: Light red/pink  
- Values = 0: White
- Values = 0.7: Light blue
- Values ≥ 1.426: Bright blue

Notes:
------
- Files use ChimeraX atom specification format (/chain:residue)
- NOT compatible with classic Chimera (use ChimeraX instead)
- Match mode is "any" to handle trimer structure (multiple chains)
- Same data applied to all 3 protomers for consistent coloring

Generated: 2025-11-21
Updated: Fixed for trimer structure with match mode "any"
