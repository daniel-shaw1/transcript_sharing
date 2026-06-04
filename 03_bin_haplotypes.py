import pandas as pd
import numpy as np
import pyranges as pr

df = pd.read_csv("cell_gene_ase.filtered.tsv", sep="\t")

gtf = pr.read_gtf("genes.gtf")

gene_pos = (
    gtf.df[gtf.df["Feature"] == "gene"]
    [["gene_name", "Chromosome", "Start", "End"]]
    .rename(columns={"gene_name": "gene"})
)

df = df.merge(gene_pos, on="gene", how="inner")

# define bins (10 Mb default)
BIN_SIZE = 10_000_000

df["bin"] = (df["Start"] // BIN_SIZE).astype(int)

df["total"] = df["B6_umi"] + df["CAST_umi"]

rows = []

for (cell, chrom, bin_id), g in df.groupby(["cell", "Chromosome", "bin"]):

    b6 = g["B6_umi"].sum()
    cast = g["CAST_umi"].sum()
    total = b6 + cast

    if total < 5:
        continue

    p = b6 / total

    rows.append([cell, chrom, bin_id, b6, cast, p])

out = pd.DataFrame(rows, columns=[
    "cell", "chrom", "bin", "B6", "CAST", "pB6"
])

out.to_csv("cell_bin_ase.tsv", sep="\t", index=False)
