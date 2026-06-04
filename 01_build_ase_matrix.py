import pysam
import pandas as pd
from collections import defaultdict

bam = pysam.AlignmentFile(
    "possorted_genome_bam.bam",
    "rb"
)

vcf = pysam.VariantFile(
    "heterozygous_snps.vcf.gz"
)

snp_lookup = {}

for rec in vcf.fetch():

    if len(rec.ref) != 1:
        continue

    if len(rec.alts) != 1:
        continue

    snp_lookup[(rec.contig, rec.pos)] = (
        rec.ref,
        rec.alts[0]
    )

counts = defaultdict(
    lambda: {"B6": set(), "CAST": set()}
)

for read in bam.fetch():

    if read.is_unmapped:
        continue

    try:
        cell = read.get_tag("CB")
        umi = read.get_tag("UB")
        gene = read.get_tag("GX")
    except KeyError:
        continue

    for qpos, rpos in read.get_aligned_pairs(matches_only=True):

        chrom = bam.get_reference_name(read.reference_id)

        key = (chrom, rpos + 1)

        if key not in snp_lookup:
            continue

        ref, alt = snp_lookup[key]

        base = read.query_sequence[qpos]

        if base == ref:
            counts[(cell, gene)]["B6"].add(umi)

        elif base == alt:
            counts[(cell, gene)]["CAST"].add(umi)

rows = []

for (cell, gene), vals in counts.items():

    rows.append(
        [
            cell,
            gene,
            len(vals["B6"]),
            len(vals["CAST"])
        ]
    )

df = pd.DataFrame(
    rows,
    columns=[
        "cell",
        "gene",
        "B6_umi",
        "CAST_umi"
    ]
)

df.to_csv(
    "cell_gene_ase.tsv",
    sep="\t",
    index=False
)
