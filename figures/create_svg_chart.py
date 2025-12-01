#!/usr/bin/env python
"""Create SVG version of the site frequency chart from build_codon_variants.ipynb cell 8."""

import tempfile
import Bio.SeqIO
import altair as alt
import dms_variants.codonvarianttable
import pandas as pd

# Disable max rows for Altair
_ = alt.data_transformers.disable_max_rows()

# Parameters from the executed notebook
gene_sequence_codon = "results/gene_sequence/codon.fasta"
gene_sequence_protein = "results/gene_sequence/protein.fasta"
site_numbering_map_csv = "data/site_numbering_map.csv"
neut_standard_barcodes = "data/neutralization_standard_barcodes.csv"

# Read neutralization standards
print(f"Reading neutralization standard barcodes from {neut_standard_barcodes}")
neut_standard = pd.read_csv(neut_standard_barcodes)

# Read gene variants
nt_variants = "results/variants/nt_variants.csv"
print(f"Reading gene variants from {nt_variants}")
gene_variants = pd.read_csv(nt_variants, na_filter=False)

# Get wildtype gene sequence
geneseq = str(Bio.SeqIO.read(gene_sequence_codon, "fasta").seq)
print(f"Read gene of length {len(geneseq)} nucleotides as primary target")

# Make a data frame with the gene variants and neutralization standard
variants_df = pd.concat(
    [
        gene_variants.assign(target="gene"),
        pd.concat(
            [
                neut_standard.rename(columns={"name": "target"}).assign(
                    library=lib,
                    substitutions="",
                    variant_call_support=gene_variants["variant_call_support"].max(),
                )
                for lib in gene_variants["library"].unique()
            ]
        ),
    ]
)

# Build the CodonVariantTables
print("Building CodonVariantTable...")
with tempfile.NamedTemporaryFile(mode="w") as f:
    variants_df.to_csv(f)
    f.flush()
    variants = dms_variants.codonvarianttable.CodonVariantTable(
        barcode_variant_file=f.name,
        geneseq=geneseq,
        allowgaps=True,
        primary_target="gene",
    )

# Read site numbering map
site_numbering_map = pd.read_csv(site_numbering_map_csv)

# Calculate site frequencies (from cell 7 in the notebook)
print("Calculating site frequencies...")
site_freqs = (
    variants.mutCounts(
        variant_type="all",
        mut_type="aa",
        samples=None,
        libraries=variants.libraries,
    )
    .query("count > 0")
    .rename(columns={"site": "sequential_site"})
    .merge(
        site_numbering_map,
        how="left",
        on="sequential_site",
        validate="many_to_one",
    )
    .assign(
        wildtype=lambda x: x["mutation"].str[0],
        mutant=lambda x: x["mutation"].str[-1],
        reference_site=lambda x: (
            x["reference_site"]
            if all(x["reference_site"] == x["reference_site"].astype(str))
            else x["reference_site"].astype("Int64")
        ),
        # remove categorical assignments to fix plotting issues
        library=lambda x: x["library"].tolist(),
    )
    .assign(mutant=lambda x: x["mutant"] + "=" + x["count"].astype(str))
    .sort_values("count", ascending=False)
    .groupby(
        ["library", "sequential_site", "reference_site", "wildtype"],
        observed=True,
        as_index=False,
    )
    .aggregate(
        count=pd.NamedAgg("count", "sum"),
        mutants=pd.NamedAgg("mutant", ", ".join),
    )
    .merge(
        variants.n_variants_df(
            samples=None, libraries=variants.libraries, primary_target_only=True
        )
        .rename(columns={"count": "n_variants"})
        .drop(columns="sample")
    )
    .assign(percent=lambda x: 100 * x["count"] / x["n_variants"])
    .drop(columns="n_variants")
)

# Create the chart (from cell 8 in the notebook)
print("Creating chart...")
zoom_brush = alt.selection_interval(
    encodings=["x"],
    mark=alt.BrushConfig(stroke="black", strokeWidth=2),
)

zoom_bar = (
    alt.Chart(site_freqs[["sequential_site", "reference_site"]].drop_duplicates())
    .mark_rect(color="lightgrey")
    .encode(
        x=alt.X(
            "reference_site",
            title=None,
            scale=alt.Scale(nice=False, zero=False),
            sort=alt.SortField("sequential_site"),
            axis=alt.Axis(labelOverlap=True),
        )
    )
    .add_params(zoom_brush)
    .properties(width=550, height=15, title="site zoom bar")
)

site_freqs_base = (
    alt.Chart()
    .encode(
        x=alt.X(
            "reference_site",
            scale=alt.Scale(nice=False, zero=False),
            sort=alt.SortField("sequential_site"),
            axis=alt.Axis(labelOverlap=True),
        ),
        y=alt.Y("percent", title="% variants with mutation"),
        tooltip=[
            (
                alt.Tooltip(c, format=".3g")
                if c in {"percent", "count", "n_variants"}
                else c
            )
            for c in site_freqs.columns
        ],
    )
    .properties(height=125, width=550)
)

site_freqs_chart = (
    alt.layer(
        site_freqs_base.mark_point(filled=True),
        site_freqs_base.mark_line(size=0.5),
        data=site_freqs,
    )
    .facet(facet=alt.Facet("library", title=None), columns=1)
    .add_params(zoom_brush)
    .transform_filter(zoom_brush)
)

site_freqs_zoom_chart = (zoom_bar & site_freqs_chart).configure_axis(grid=False)

# Save as SVG
output_file = "site_freqs_chart.svg"
print(f"Saving chart to {output_file}...")
site_freqs_zoom_chart.save(output_file)
print(f"Successfully saved SVG to {output_file}")
