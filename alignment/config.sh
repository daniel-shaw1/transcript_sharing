# Shared settings for the alignment workflow. Source this file; do not run it.

# Where everything is written (needs ~400 GB for all E-MTAB-11602 FASTQs + index)
WORKDIR="${WORKDIR:-$PWD/work}"

THREADS="${THREADS:-16}"

# ArrayExpress accession for Panten et al., Nat Commun 15, 1272 (2024)
# (E-MTAB-6934 holds the earlier B6-only samples, not needed for F1 ASE)
STUDY="E-MTAB-11602"

# mm10 / GRCm38, Ensembl 102 -- same coordinates as the Sanger MGP REL-1505 SNP calls
GENOME_FA_URL="https://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz"
GTF_URL="https://ftp.ensembl.org/pub/release-102/gtf/mus_musculus/Mus_musculus.GRCm38.102.gtf.gz"

# Sanger Mouse Genomes Project SNPs (mm10). Ref allele = C57BL/6J, alt = CAST/EiJ
MGP_VCF_URL="https://ftp.ebi.ac.uk/pub/databases/mousegenomes/REL-1505-SNPs_Indels/mgp.v5.merged.snps_all.dbSNP142.vcf.gz"
CAST_SAMPLE="CAST_EiJ"

# 10x barcode whitelists (v3 and v2 chemistry)
WL_V3_URL="https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/3M-february-2018.txt.gz"
WL_V2_URL="https://raw.githubusercontent.com/10XGenomics/cellranger/master/lib/python/cellranger/barcodes/737K-august-2016.txt"

REF_DIR="$WORKDIR/ref"
FASTQ_DIR="$WORKDIR/fastq"
ALIGN_DIR="$WORKDIR/aligned"
