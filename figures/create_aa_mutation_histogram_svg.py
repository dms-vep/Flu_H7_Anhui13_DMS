#!/usr/bin/env python
"""Create SVG version of the amino acid mutation histogram from build_codon_variants.ipynb."""

import tempfile
import Bio.SeqIO
import dms_variants.codonvarianttable
import pandas as pd

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

# Create the histogram plot
print("Creating amino acid mutation histogram...")
p = variants.plotNumMutsHistogram(
    mut_type="aa",
    max_muts=8,
    widthscale=1.3,
    heightscale=1,
    libraries=variants.libraries,
    samples=None,
)
fig = p.draw(show=False)

# Save as SVG
output_file = "aa_mutation_histogram.svg"
print(f"Saving figure to {output_file}...")
fig.savefig(output_file, format='svg', bbox_inches='tight')
print(f"Successfully saved SVG to {output_file}")
