---
aside: false
---

# Cell Entry

We measured how mutations to the H7 HA affected its ability to mediate entry into 293 cells expressing only a2-3-linked sialic acids, only a2-6 linked sialic acids, or an equal mix of the two cells.
Below are interactive plots showing the results and links to CSVs with the data.

Note that H3 numbering is used for all plots (see [here](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4229193/)).

## Plot showing entry in all three types of cells
Below is a plot showing the effects of mutations to the HA ectodomain on entry into all three cells.
A standalone version of this plot is available [here](htmls/entry_all_cells_overlaid.html){target="_self"}.
The data in this plot is filtered to remove low quality measurements.

<Figure caption="Effects of mutations on H7 HA entry in 293 cells expressing a2-3 sialic acids, a2-6 sialic acids, or a mix of both types of cells">
    <Altair :showShadow="true" :spec-url="'htmls/entry_all_cells_overlaid.html'"></Altair>
</Figure>

## Scatter plot comparing effects of mutations in all three types of cells
Below is a plot comparing the average effect of all mutations at each site in HA on entry into all three cells; click [here](htmls/cell_entry_scatter.html){target="_self"} for a standalone version of this plot.

<Figure caption="Scatter plot comparing average effects of mutations at each site in all three cell types">
    <Altair :showShadow="true" :spec-url="'htmls/cell_entry_scatter.html'"></Altair>
</Figure>

## Separate plots for entry in each type of cell
The links below are for interactive plots for H7 HA-mediated entry into each cell type plotted separately, along with detailed QC metrics that can be seen by mousing over points (eg, per-replicate measurements) and adjusting the sliders:

 - [entry in 293 cells expressing only a2-3 linked sialic acids](htmls/293_2-3_entry_func_effects.html){target="_self"}
 - [entry in 293 cells expressing only a2-6 linked sialic acids](htmls/293_2-6_entry_func_effects.html){target="_self"}
 - [entry in an equal mix of the a2-3 and a2-6 293 cells](htmls/293_mix_entry_func_effects.html){target="_self"}

## Numerical data
Click [here](https://github.com/dms-vep/Flu_H7_Anhui13_DMS/blob/main/results/summaries/entry_all_cells.csv) for a CSV giving the effects of mutations on entry into each cell.
These are measurements after filtering away low-quality measurements.

See the [Appendix](appendix.html){target="_self"} for links to additional files with values prior to QC and with additional QC metrics.

