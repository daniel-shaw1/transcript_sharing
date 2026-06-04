import pandas as pd
import numpy as np

df = pd.read_csv("cell_gene_ase.tsv", sep="\t")

df["total"] = df["B6_umi"] + df["CAST_umi"]

# basic gene-level filter
df = df[df["total"] >= 2]

# cell-level QC
cell_qc = df.groupby("cell")["total"].sum().reset_index()
cell_qc.columns = ["cell", "cell_umi"]

good_cells = cell_qc[cell_qc["cell_umi"] >= 1000]["cell"]

df = df[df["cell"].isin(good_cells)]

# gene-level QC
gene_qc = df.groupby("gene")["total"].sum().reset_index()
gene_qc.columns = ["gene", "gene_umi"]

good_genes = gene_qc[gene_qc["gene_umi"] >= 50]["gene"]

df = df[df["gene"].isin(good_genes)]

# spermatid enrichment step (placeholder)
# assume you already have annotation table
anno = pd.read_csv("cell_annotations.tsv", sep="\t")
spermatid_cells = anno[anno["celltype"].isin([
    "round_spermatid",
    "elongating_spermatid"
])]["cell"]

df = df[df["cell"].isin(spermatid_cells)]

df.to_csv("cell_gene_ase.filtered.tsv", sep="\t", index=False)
