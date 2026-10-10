# Alignment of Panten et al. (2024) 10x data to mm10

Data: Panten et al., *The dynamic genetic determinants of increased transcriptional
divergence in spermatids*, Nat Commun 15, 1272 (2024),
[doi:10.1038/s41467-024-45133-1](https://doi.org/10.1038/s41467-024-45133-1).
Raw 10x 3′ FASTQs are in ArrayExpress **E-MTAB-11602** (B6, CAST, and B6×CAST F1 testis).
The authors' own code is at [PMBio/ase_spermatogenesis](https://github.com/PMBio/ase_spermatogenesis);
they used Cell Ranger defaults plus WASP-based allele counting.

This workflow produces `possorted_genome_bam.bam` and `heterozygous_snps.vcf.gz`,
the two inputs `01_build_ase_matrix.py` reads.

| Choice | Value |
|---|---|
| Genome | GRCm38 / mm10, Ensembl 102 primary assembly |
| Annotation | Ensembl 102 GTF, Cell Ranger biotype filter |
| SNPs | Sanger MGP REL-1505, CAST/EiJ hom-alt PASS SNPs (ref = B6, alt = CAST) |
| Aligner | STARsolo, Cell Ranger-equivalent settings, chemistry from R1 length |
| Mapping bias | STAR WASP mode; reads failing WASP (`vW != 1`) removed |

## Requirements

Run this on an HPC node, not a laptop or the cloud session:
- ≥ 32 GB RAM (STAR mouse index), 16 cores recommended
- ~400 GB disk for all libraries; much less if only the F1 libraries are kept
- network access to `www.ebi.ac.uk`, `ftp.ebi.ac.uk`, `ftp.sra.ebi.ac.uk`,
  `ftp.ensembl.org`, `raw.githubusercontent.com`

## Running

```bash
conda env create -f alignment/environment.yml
conda activate transcript_sharing_align
export WORKDIR=/scratch/$USER/transcript_sharing THREADS=16

bash alignment/00_fetch_data.sh            # writes $WORKDIR/fastq/runs.tsv
# edit runs.tsv: keep the F1 B6xCAST libraries only
bash alignment/00_fetch_data.sh download
bash alignment/01_build_reference.sh
bash alignment/02_align_starsolo.sh all
```

Then, per F1 sample, from that sample's output directory:

```bash
cd $WORKDIR/aligned/<sample>
ln -s $WORKDIR/ref/heterozygous_snps.vcf.gz* .
python /path/to/transcript_sharing/01_build_ase_matrix.py
```
