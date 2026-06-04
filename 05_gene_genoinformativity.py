import pandas as pd
import numpy as np
from scipy.optimize import minimize

ase = pd.read_csv("cell_gene_ase.filtered.tsv", sep="\t")
hap = pd.read_csv("cell_bin_haplotypes.tsv", sep="\t")

# map gene → bin (reuse GTF mapping logic)
gene_map = pd.read_csv("gene_to_bin.tsv", sep="\t")

df = ase.merge(gene_map, on="gene")
df = df.merge(
    hap[["cell", "chrom", "bin", "H"]],
    on=["cell", "chrom", "bin"],
    how="inner"
)

df["total"] = df["B6_umi"] + df["CAST_umi"]

df = df[df["total"] >= 2]


def neg_log_likelihood(g, data):

    g = g[0]

    g = np.clip(g, 0, 1)

    ll = 0

    for row in data.itertuples():

        H = row.H
        total = row.total
        b6 = row.B6_umi

        p_h = 0.5 + 0.5 * g * (2 * H - 1)

        p_h = np.clip(p_h, 1e-3, 1 - 1e-3)

        ll += b6 * np.log(p_h) + (total - b6) * np.log(1 - p_h)

    return -ll


results = []

for gene, gdf in df.groupby("gene"):

    if len(gdf) < 20:
        continue

    res = minimize(
        neg_log_likelihood,
        x0=[0.1],
        args=(gdf,),
        bounds=[(0, 1)]
    )

    results.append([
        gene,
        res.x[0],
        len(gdf)
    ])

out = pd.DataFrame(results, columns=[
    "gene", "g", "n_obs"
])

out.sort_values("g", ascending=False).to_csv(
    "gene_genoinformativity.tsv",
    sep="\t",
    index=False
)
