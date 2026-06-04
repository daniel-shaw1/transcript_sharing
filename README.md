# transcript_sharing

# Genoinformativity from 10x Single-Cell Spermatids

This repository implements a pipeline to estimate **gene-level genoinformativity (transcript sharing sensitivity to haplotype)** from 10x Genomics single-cell RNA-seq data in F1 hybrid spermatogenesis datasets (e.g., B6 × CAST).

The method adapts ideas from Bhutani et al. (Science, 2021) to sparse 10x data by leveraging:

- large cell numbers
- allele-specific UMIs
- recombination-aware haplotype inference
- bin-level smoothing instead of per-SNP inference

---

## Biological goal

We estimate, for each gene:

g_j ∈ [0,1]

where:

- g = 0 → complete transcript sharing between spermatids  
- g = 1 → no sharing (fully genotype-informative transcripts)  
- intermediate values → partial retention / limited sharing  

This captures how strongly a gene’s expression reflects its underlying haplotype.

---

## Key idea

Instead of modeling transcript sharing directly, we measure:

> How predictable is local haplotype from allele-specific expression?

Recombination in spermatids provides natural “labels” for local genotype inference.

---

## Pipeline overview

FASTQ
  ↓
Cell Ranger / STARsolo
  ↓
Allele-specific UMI extraction (VCF-based)
  ↓
QC + spermatid filtering
  ↓
Gene × cell ASE matrix
  ↓
Chromosomal bin aggregation
  ↓
Haplotype inference (HMM over bins)
  ↓
Gene-level likelihood model
  ↓
Genoinformativity (g) estimation

---

## Data requirements

### Required inputs

- 10x FASTQs (B6 × CAST F1 or similar cross)
- reference genome (mm10/mm39)
- parental SNP VCF (CAST/B6 or equivalent)
- gene annotation GTF

### Assumptions

- diploid F1 hybrid system
- sufficient heterozygous SNP density
- standard 10x 3′ scRNA-seq chemistry

---

## Repository structure

gim10x/
├── scripts/
│   ├── 02_filter_cells.py
│   ├── 03_bin_haplotypes.py
│   ├── 04_hmm_haplotypes.py
│   ├── 05_gene_genoinformativity.py
│   └── utils.py
├── data/
├── results/
└── README.txt

---

## Step 1 — QC and filtering

Script: 02_filter_cells.py

Filters:
- low UMI cells (<1000 UMIs)
- low-information genes (<50 total UMIs)
- non-spermatid cell types
- genes with insufficient ASE signal

Output:
cell_gene_ase.filtered.tsv

---

## Step 2 — Gene-level ASE aggregation

Allele-specific UMIs are aggregated per gene:

cell | gene | B6_umi | CAST_umi

Only SNP-informative UMIs are counted.

---

## Step 3 — Chromosomal binning

Script: 03_bin_haplotypes.py

Genes are mapped to genomic bins (default: 10 Mb). ASE is aggregated per:

- cell
- chromosome
- bin

Output:
cell_bin_ase.tsv

This reduces sparsity and enables recombination-aware inference.

---

## Step 4 — Haplotype inference (recombination-aware)

Script: 04_hmm_haplotypes.py

Each cell is modeled as a mosaic of parental haplotypes.

Model:
- States: B6, CAST
- Observations: bin-level allelic fractions
- Transitions: recombination events
- Output: haplotype blocks per cell

Output:
cell_bin_haplotypes.tsv

---

## Step 5 — Genoinformativity estimation

Script: 05_gene_genoinformativity.py

For each gene:

P(B6 reads) =
0.5 + (g/2)(2H - 1)

where:
- H = inferred local haplotype
- g = genoinformativity parameter

Fitting:
- maximum likelihood per gene
- g constrained to [0,1]

Output:
gene_genoinformativity.tsv

Columns:
- gene
- g
- n_obs

---

## Interpretation

g value → meaning

~0.0 → full transcript sharing  
~0.3 → weak genotype dependence  
~0.6 → strong partial retention  
~1.0 → near-complete genotype encoding  

---

## Key assumptions

- ASE signal is not dominated by mapping bias
- cis-eQTL effects are not fully confounded with sharing
- haplotype inference is sufficiently accurate at bin level
- recombination provides sufficient mixing across cells
- F1 heterozygosity is known and provided via VCF

---

## Limitations

- 10x sparsity limits per-gene precision
- cis-regulatory variation can inflate g estimates
- haplotype inference is approximate (bin-level)
- low-expression genes are unreliable
- assumes F1 hybrid system

---

## Recommended improvements

- Beta-binomial emission model instead of binomial
- Forward-backward HMM instead of Viterbi
- Mutual information-based genoinformativity
- Hierarchical shrinkage across gene classes
- Pseudotime-dependent g_j(t)
- Explicit cis-eQTL vs sharing decomposition

---

## Quick start

# 1. filter ASE matrix
python scripts/02_filter_cells.py

# 2. bin haplotypes
python scripts/03_bin_haplotypes.py

# 3. infer haplotypes
python scripts/04_hmm_haplotypes.py

# 4. compute genoinformativity
python scripts/05_gene_genoinformativity.py

---

## Citation

If you use this pipeline, please consider citing:

- Bhutani et al., Science (2021)
- 10x Genomics single-cell RNA-seq methodology papers
- relevant F1 hybrid ASE literature

---

## Notes

