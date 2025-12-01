#!/usr/bin/env python
"""Create SVG version of the mutation distribution chart from build_codon_variants.ipynb."""

import tempfile
import Bio.SeqIO
import altair as alt
import dms_variants.codonvarianttable
import pandas as pd

# Disable max rows for Altair
_ = alt.data_transformers.disable_max_rows()

# Parameters from the executed notebook
gene_sequence_codon = "results/gene_sequence/codon.fasta"
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

# Create the mutation distribution data (from cell 8)
print("Calculating mutation distribution...")
num_muts = variants.numCodonMutsByType(
    variant_type="all",
    samples=None,
    libraries=variants.libraries,
).assign(
    # remove categorical assignments to fix plotting issues
    mutation_type=lambda x: x["mutation_type"].tolist(),
    library=lambda x: x["library"].tolist(),
)

mut_type_order = (
    num_muts.sort_values("number", ascending=False)["mutation_type"].unique().tolist()
)

# Create the chart
print("Creating chart...")
num_muts_chart = (
    alt.Chart(num_muts)
    .encode(
        x=alt.X("mutation_type", title=None, scale=alt.Scale(domain=mut_type_order)),
        y=alt.Y("number", title="codon mutations per variant"),
        color=alt.Color("mutation_type", legend=None),
        facet=alt.Facet("library", title=None, columns=4),
        tooltip=["library", "mutation_type", alt.Tooltip("number", format=".2f")],
    )
    .mark_bar()
    .properties(width=alt.Step(25), height=150)
)

# Save as SVG
output_file = "mutation_distribution_chart.svg"
print(f"Saving chart to {output_file}...")
num_muts_chart.save(output_file)
print(f"Successfully saved SVG to {output_file}")
