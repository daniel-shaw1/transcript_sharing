import pandas as pd
import numpy as np

df = pd.read_csv("cell_bin_ase.tsv", sep="\t")

TRANSITION = 1e-3  # recombination prior

def viterbi(chrom_df):

    chrom_df = chrom_df.sort_values("bin")

    n = len(chrom_df)

    states = ["B6", "CAST"]

    logp = np.zeros((n, 2))
    path = np.zeros((n, 2), dtype=int)

    for i, row in enumerate(chrom_df.itertuples()):

        p = row.pB6
        p = min(max(p, 1e-3), 1 - 1e-3)

        # likelihoods
        loglik_b6 = np.log(p)
        loglik_cast = np.log(1 - p)

        if i == 0:
            logp[i] = [loglik_b6, loglik_cast]
            continue

        for s in range(2):

            stay = logp[i-1, s] + np.log(1 - TRANSITION)
            switch = logp[i-1, 1-s] + np.log(TRANSITION)

            emission = loglik_b6 if s == 0 else loglik_cast

            logp[i, s] = max(stay, switch) + emission

    # backtrack (simple version)
    states_out = np.zeros(n)

    states_out[-1] = np.argmax(logp[-1])

    for i in range(n-2, -1, -1):
        states_out[i] = np.argmax(logp[i])

    chrom_df["H"] = states_out

    return chrom_df


out = []

for (cell, chrom), g in df.groupby(["cell", "chrom"]):

    res = viterbi(g)

    out.append(res)

out = pd.concat(out)

out.to_csv("cell_bin_haplotypes.tsv", sep="\t", index=False)
