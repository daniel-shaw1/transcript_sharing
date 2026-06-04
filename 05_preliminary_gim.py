import pandas as pd
import numpy as np

df = pd.read_csv(
    "cell_gene_ase.tsv",
    sep="\t"
)

df["total"] = (
    df["B6_umi"] +
    df["CAST_umi"]
)

df = df[df["total"] >= 2]

df["AI"] = (
    df["B6_umi"] /
    df["total"]
)

stats = []

for gene, g in df.groupby("gene"):

    if len(g) < 30:
        continue

    var = g["AI"].var()

    gi = min(
        1.0,
        4 * var
    )

    stats.append(
        [gene, len(g), gi]
    )

out = pd.DataFrame(
    stats,
    columns=[
        "gene",
        "n_cells",
        "GI"
    ]
)

out.sort_values(
    "GI",
    ascending=False
).to_csv(
    "preliminary_gim.tsv",
    sep="\t",
    index=False
)
