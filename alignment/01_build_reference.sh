#!/usr/bin/env bash
# Build the mm10 STAR index and the B6/CAST SNP VCF used for WASP filtering
# and allele calling (01_build_ase_matrix.py).
#
# Needs ~32 GB RAM for genomeGenerate.
set -euo pipefail
source "$(dirname "$0")/config.sh"
mkdir -p "$REF_DIR"
cd "$REF_DIR"

[[ -s genome.fa ]] || curl -fSL "$GENOME_FA_URL" | gunzip -c > genome.fa
[[ -s genes.gtf ]] || curl -fSL "$GTF_URL" | gunzip -c > genes.gtf
[[ -s genome.fa.fai ]] || samtools faidx genome.fa

# Keep genes with biotypes Cell Ranger uses, so counts match the paper's processing
if [[ ! -s genes.filtered.gtf ]]; then
    grep -E '^#|gene_biotype "(protein_coding|lncRNA|IG_[A-Z]*_gene|TR_[A-Z]*_gene)"' genes.gtf \
        > genes.filtered.gtf
fi

# CAST/EiJ homozygous-alt PASS SNPs. Ref allele = B6, alt allele = CAST, which is the
# convention 01_build_ase_matrix.py assumes. Genotype is rewritten to 0|1 because the
# F1 is heterozygous at every such site and STAR's WASP mode only uses het sites.
if [[ ! -s heterozygous_snps.vcf.gz ]]; then
    [[ -s mgp.snps.vcf.gz ]] || curl -fSL -o mgp.snps.vcf.gz "$MGP_VCF_URL"
    [[ -s mgp.snps.vcf.gz.tbi ]] || curl -fSL -o mgp.snps.vcf.gz.tbi "$MGP_VCF_URL.tbi"
    bcftools view -s "$CAST_SAMPLE" -v snps -m2 -M2 -f PASS -i 'GT="1/1"' mgp.snps.vcf.gz -Ou \
        | bcftools annotate -x INFO,^FORMAT/GT -Ov \
        | awk 'BEGIN{OFS="\t"} /^#/ {print; next} {$10="0|1"; print}' \
        | bcftools view -Oz -o heterozygous_snps.vcf.gz
    bcftools index -t heterozygous_snps.vcf.gz
fi
# Plain VCF for STAR, which cannot read gzip
[[ -s heterozygous_snps.vcf ]] || bcftools view heterozygous_snps.vcf.gz -Ov -o heterozygous_snps.vcf

# 10x whitelists
[[ -s whitelist_v3.txt ]] || curl -fsSL "$WL_V3_URL" | gunzip -c > whitelist_v3.txt
[[ -s whitelist_v2.txt ]] || curl -fsSL -o whitelist_v2.txt "$WL_V2_URL"

if [[ ! -s star_index/SA ]]; then
    mkdir -p star_index
    STAR --runMode genomeGenerate \
        --runThreadN "$THREADS" \
        --genomeDir star_index \
        --genomeFastaFiles genome.fa \
        --sjdbGTFfile genes.filtered.gtf \
        --sjdbOverhang 100
fi

echo "Reference ready in $REF_DIR"
