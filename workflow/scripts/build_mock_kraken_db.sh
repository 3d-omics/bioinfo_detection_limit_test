#!/usr/bin/env bash
set -euo pipefail

# Assume kraken2 is installed

# Build db


# Has to be a fasta
kraken2-build \
    --no-masking \
    --add-to-library mags_sub_kraken.fa \
    --db resources/kraken_mock \
    --threads 4

kraken2-build \
    --download-taxonomy \
    --db resources/kraken_mock

kraken2-build \
    --build \
    --db resources/kraken_mock

kraken2-build \
    --clean \
    --db resources/kraken_mock

# Test

kraken2 \
    --db resources/kraken_mock \
    --threads 4 \
    --report resources/kraken_mock.report \
    --output resources/kraken_mock.out \
    --gzip-compressed \
    resources/reference/mags_sub_kraken.fa.gz


# Run

kraken2 \
    --db resources/kraken_mock \
    --threads 8 \
    --report results/kraken/sample1.lib1.report \
    --output resources/kraken/sample1.lib1.out \
    --gzip-compressed \
    --paired \
    results/fastp/sample1.lib1_1.fq.gz \
    results/fastp/sample1.lib1_2.fq.gz


kraken2 \
    --db resources/kraken_mock \
    --threads 8 \
    --report results/kraken/sample2.lib1.report \
    --output resources/kraken/sample2.lib1.out \
    --gzip-compressed \
    --paired \
    results/fastp/sample2.lib1_1.fq.gz \
    results/fastp/sample2.lib1_2.fq.gz
