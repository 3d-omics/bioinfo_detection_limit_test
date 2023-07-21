#!/usr/bin/env bash
set -euo pipefail

>&2 echo -e "Make sure multiqc, coverm, r-tidyverse, r-argparse and r-nonpareil are in your \$PATH"
>&2 echo -e "mamba create -n failed_reports coverm multiqc r-tidyverse r-argparse r-nonpareil"
>&2 echo -e "conda activate failed_reports"


# Get kraken2 stats
multiqc \
    --title kraken2 \
    --force \
    --filename kraken2 \
    --outdir reports/by_step \
    results/kraken2/*.report \
2> results/by_step/kraken2.log 1>&2

# Get nonpareil
Rscript --no-init-file workflow/scripts/aggregate_nonpareil.R \
    --input-folder results/stats/nonpareil/ \
    --output-file results/stats/nonpareil.tsv \
2> resiults/stats/nonpareil.log 1>&2


# Get singlem
singlem condense \
    --input-archive-otu-tables results/stats/singlem/*.archive.json \
    --taxonomic-profile results/stats/singlem.tsv \
    --metapackage results/stats/singlem/data/S3.2.0.GTDB_r214.metapackage_20230428.smpkg.zb \
2> results/stats/singlem.log 1>&2


# Get coverm
coverm genome \
    --bam-files results/stats/coverm/*.bam \
    --methods count \
    --separator ^ \
    --threads 8 \
    --min-covered-fraction 0.0 \
> results/stats/coverm_genome.tsv \
2> results/stats/coverm_genome.log


coverm contig \
    --threads 8 \
    --bam-files results/stats/coverm/*.bam \
    --methods count \
    --proper-pairs-only \
> results/stats/coverm_contig.tsv \
2> results/stats/coverm_contig.log
