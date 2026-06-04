import numpy as np
import pandas as pd
from scipy.optimize import minimize

# -----------------------
# LOAD DATA
# -----------------------

ase = pd.read_csv("cell_gene_ase.tsv", sep="\t")
hap = pd.read_csv("cell_bin_haplotypes.tsv", sep="\t")

ase["total"] = ase["B6_umi"] + ase["CAST_umi"]
ase = ase[ase["total"] > 0]


# -----------------------
# BUILD SIMPLE CELL GRAPH
# (must match simulator structure)
# -----------------------

def adjacency_matrix(n, window=5):
    A = np.zeros((n, n))
    for i in range(n):
        for j in range(max(0, i-window), min(n, i+window+1)):
            if i != j:
                A[i, j] = 1
    A = A / A.sum(axis=1, keepdims=True)
    return A


cells = sorted(ase["cell"].unique())
cell_index = {c:i for i,c in enumerate(cells)}

A = adjacency_matrix(len(cells), window=5)


# -----------------------
# MAP HAPLOTYPE TO MATRIX
# -----------------------

hap["cell_i"] = hap["cell"].map(cell_index)
hap = hap.dropna()
hap["cell_i"] = hap["cell_i"].astype(int)


# pivot: cell x gene haplotype
H = hap.pivot(index="cell_i", columns="bin", values="H").fillna(0).values


# -----------------------
# DIFFUSION OPERATOR
# -----------------------

def diffuse(H, A, alpha):

    return (1 - alpha) * H + alpha * (A @ H)


# -----------------------
# NEG LOG LIKELIHOOD
# -----------------------

def nll(params, gdf, H_eff):

    alpha = params[0]
    g = params[1]

    alpha = np.clip(alpha, 0, 1)
    g = np.clip(g, 0, 1)

    ll = 0.0

    for row in gdf.itertuples():

        c = cell_index[row.cell]
        total = row.total
        b6 = row.B6_umi

        h = H_eff[c]

        p = 0.5 + 0.5 * g * (2*h - 1)
        p = np.clip(p, 1e-6, 1-1e-6)

        ll += b6*np.log(p) + (total-b6)*np.log(1-p)

    return -ll


# -----------------------
# FIT PER GENE
# -----------------------

results = []

genes = ase["gene"].unique()

for gene in genes:

    gdf = ase[ase["gene"] == gene]

    if len(gdf) < 20:
        continue

    # initial guess
    params0 = np.array([0.3, 0.3])  # alpha, g

    def objective(params):
        H_eff = diffuse(H[:,0], A, params[0])  # simplified per-gene proxy
        return nll(params, gdf, H_eff)

    res = minimize(
        objective,
        params0,
        bounds=[(0,1),(0,1)],
        method="L-BFGS-B"
    )

    results.append([gene, res.x[0], res.x[1], len(gdf)])


out = pd.DataFrame(results, columns=[
    "gene", "alpha_hat", "g_hat", "n_obs"
])

out.to_csv("diffusion_g_results.tsv", sep="\t", index=False)

print("Done.")
