#!/usr/bin/env bash
# Download the E-MTAB-11602 10x FASTQs (Panten et al. 2024).
#
#   bash alignment/00_fetch_data.sh            # list runs only -> runs.tsv
#   bash alignment/00_fetch_data.sh download   # also download FASTQs in runs.tsv
#
# Edit runs.tsv between the two calls to keep only the F1 (B6 x CAST / CAST x B6)
# libraries; the parental B6 / CAST libraries carry no allele-specific signal.
set -euo pipefail
source "$(dirname "$0")/config.sh"
mkdir -p "$FASTQ_DIR"

ENA="https://www.ebi.ac.uk/ena/portal/api"
RUNS="$FASTQ_DIR/runs.tsv"

if [[ ! -s "$RUNS" ]]; then
    # ArrayExpress studies are mirrored to ENA with the E-MTAB id as the study alias
    PROJECT=$(curl -fsS "$ENA/search?result=study&query=study_alias=%22$STUDY%22&fields=study_accession&format=tsv" \
        | awk 'NR==2 {print $1}')
    if [[ -z "$PROJECT" ]]; then
        echo "Could not resolve $STUDY to an ENA project; set PROJECT=PRJEBxxxx by hand" >&2
        exit 1
    fi
    echo "$STUDY -> $PROJECT"
    curl -fsS "$ENA/filereport?accession=$PROJECT&result=read_run&format=tsv&fields=run_accession,sample_alias,sample_title,library_name,submitted_ftp,fastq_ftp" \
        > "$RUNS"
    echo "Wrote $(($(wc -l < "$RUNS") - 1)) runs to $RUNS"
    echo "Keep only the F1 libraries in $RUNS, then rerun with 'download'."
fi

[[ "${1:-}" == "download" ]] || exit 0

# Prefer the submitter's original files (10x I1/R1/R2 FASTQs); ENA's own
# fastq_ftp sometimes holds only R2 for 10x submissions.
tail -n +2 "$RUNS" | while IFS=$'\t' read -r run sample_alias _ _ submitted fastq; do
    files="${submitted:-$fastq}"
    dest="$FASTQ_DIR/$sample_alias"
    mkdir -p "$dest"
    for f in ${files//;/ }; do
        out="$dest/$(basename "$f")"
        [[ -s "$out" ]] && continue
        echo "[$run] $f"
        curl -fSL --retry 5 -C - -o "$out.part" "https://$f" && mv "$out.part" "$out"
    done
done

# Processed matrices + the authors' allele-specific counts (useful to cross-check)
mkdir -p "$WORKDIR/processed"
curl -fsS "https://ftp.ebi.ac.uk/biostudies/fire/E-MTAB-/${STUDY: -3}/$STUDY/Files/" \
    -o "$WORKDIR/processed/index.html" || echo "BioStudies listing unavailable; skipping processed files"
