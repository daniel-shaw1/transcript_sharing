#!/usr/bin/env bash
# Align one 10x library to mm10 with STARsolo, using Cell Ranger-equivalent
# settings plus WASP filtering against the B6/CAST SNPs.
#
#   bash alignment/02_align_starsolo.sh <sample_dir>     # e.g. work/fastq/Sample5
#   bash alignment/02_align_starsolo.sh all              # every dir under $FASTQ_DIR
#
# Output per sample (in $ALIGN_DIR/<sample>/):
#   possorted_genome_bam.bam(.bai)   coord-sorted, WASP-passing reads, CB/UB/GX/GN tags
#   Solo.out/Gene/filtered/          cell x gene matrix (Cell Ranger format)
# The BAM is what 01_build_ase_matrix.py expects.
set -euo pipefail
source "$(dirname "$0")/config.sh"

align_sample() {
    local sample_dir="$1"
    local sample
    sample=$(basename "$sample_dir")
    local out="$ALIGN_DIR/$sample"

    if [[ -s "$out/possorted_genome_bam.bam.bai" ]]; then
        echo "[$sample] already aligned, skipping"
        return
    fi
    mkdir -p "$out"

    # Cell Ranger names (*_R1_001.fastq.gz) or ENA names (*_1.fastq.gz)
    local r1 r2
    r1=$(ls "$sample_dir"/*_R1_*.fastq.gz "$sample_dir"/*_1.fastq.gz 2>/dev/null | sort | paste -sd, || true)
    r2=$(ls "$sample_dir"/*_R2_*.fastq.gz "$sample_dir"/*_2.fastq.gz 2>/dev/null | sort | paste -sd, || true)
    if [[ -z "$r1" || -z "$r2" ]]; then
        echo "[$sample] no R1/R2 FASTQs found in $sample_dir" >&2
        return 1
    fi

    # Chemistry from barcode read length: 26 bp = v2 (10 bp UMI), >=28 bp = v3 (12 bp UMI)
    local r1_len umi_len whitelist
    r1_len=$(zcat "${r1%%,*}" | head -2 | tail -1 | tr -d '\n' | wc -c)
    if (( r1_len >= 28 )); then
        umi_len=12; whitelist="$REF_DIR/whitelist_v3.txt"
    else
        umi_len=10; whitelist="$REF_DIR/whitelist_v2.txt"
    fi
    echo "[$sample] R1 length $r1_len -> UMI length $umi_len"

    STAR --runThreadN "$THREADS" \
        --genomeDir "$REF_DIR/star_index" \
        --readFilesIn "$r2" "$r1" \
        --readFilesCommand zcat \
        --outFileNamePrefix "$out/" \
        --soloType CB_UMI_Simple \
        --soloCBwhitelist "$whitelist" \
        --soloCBstart 1 --soloCBlen 16 \
        --soloUMIstart 17 --soloUMIlen "$umi_len" \
        --soloBarcodeReadLength 0 \
        --soloCBmatchWLtype 1MM_multi_Nbase_pseudocounts \
        --soloUMIfiltering MultiGeneUMI_CR \
        --soloUMIdedup 1MM_CR \
        --soloCellFilter EmptyDrops_CR \
        --soloFeatures Gene GeneFull \
        --clipAdapterType CellRanger4 \
        --outFilterScoreMin 30 \
        --waspOutputMode SAMtag \
        --varVCFfile "$REF_DIR/heterozygous_snps.vcf" \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM vA vG vW \
        --limitBAMsortRAM 30000000000

    # Drop reads that overlap a SNP but fail WASP (vW != 1); keep reads without SNPs
    samtools view -@ "$THREADS" -b \
        -e '!([vW]) || [vW]==1' \
        -o "$out/possorted_genome_bam.bam" \
        "$out/Aligned.sortedByCoord.out.bam"
    samtools index -@ "$THREADS" "$out/possorted_genome_bam.bam"
    rm "$out/Aligned.sortedByCoord.out.bam"

    echo "[$sample] done: $out/possorted_genome_bam.bam"
}

if [[ "${1:-}" == "all" ]]; then
    for d in "$FASTQ_DIR"/*/; do
        align_sample "${d%/}"
    done
elif [[ -n "${1:-}" ]]; then
    align_sample "$1"
else
    echo "usage: $0 <sample_fastq_dir> | all" >&2
    exit 1
fi
