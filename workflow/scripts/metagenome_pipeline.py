#!/usr/bin/env python

"""
===============================================================================
Python script name metagenome.pipeline.py It's essential part of the pipeline
and executed from metagenome.pipeline.sh bash script
Script written by Alexander Kukalev
Ana Pombo lab, Berlin Institute of Medical systems Biology,
Max-Delbrueck Center for Molecular Medizine, Berlin, Germany
Date: August 25, 2023

This pipeline was developed as part of 3DOmics project. It takes any number of
fastq raw sequencing files and compute the reads and coverage for different
bacterial species, as well as human and chicken contaminations to determine
quality of DNA extraction and quality metrics.

It requires the following software package to be installed:
- bedtools
- samtools
- bowtie2
- pandas

This can be done with the following conda command:

conda create \
    -n 3domics_pombo \
    -c bioconda \
    bedtools samtools bowtie2 pandas

Usage example:
$./metagenome.pipeline.sh input_folder name_of_the_dataset

===============================================================================
"""

import sys

import pandas as pd


def remove_path_from_columns_names(data_frame):
    """remove folder names and extension from the column names"""
    list_of_columns = data_frame.columns.tolist()
    simplified_columns = []
    for file_name in list_of_columns:
        file_name = file_name.split("/")[-1]
        file_name = file_name.split(".")[0]
        if file_name[0].isdigit():  # Add letter to the file name is starts from digit
            file_name = "f_" + file_name
        file_name = file_name.replace("-", "_")
        simplified_columns.append(file_name)
    data_frame.columns = simplified_columns
    return data_frame


if __name__ == "__main__":
    READS_TABLE_PATH = sys.argv[1]
    COVERAGE_TABLE_PATH = sys.argv[2]
    DATASET_NAME = sys.argv[3]

    reads_df = pd.read_csv(READS_TABLE_PATH, sep="\t")
    coverage_df = pd.read_csv(COVERAGE_TABLE_PATH, sep="\t")

    # Simplify column names
    reads_counts_df = remove_path_from_columns_names(reads_df)
    coverage_counts_df = remove_path_from_columns_names(coverage_df)

    # Remove coordinates columns
    reads_counts_df = reads_counts_df.drop(columns=["start", "stop"])
    coverage_counts_df = coverage_counts_df.drop(columns=["start", "stop"])

    # Aggregate species
    reads_counts_df.loc[:, "specie"] = reads_counts_df["chrom"].str.split("_").str[0]
    coverage_counts_df.loc[:, "specie"] = (
        coverage_counts_df["chrom"].str.split("_").str[0]
    )
    reads_counts_sum = reads_counts_df.groupby(["specie"]).sum()
    coverage_counts_sum = coverage_counts_df.groupby(["specie"]).sum()

    # Calculate percentage for reads
    total_reads = reads_counts_sum.sum()
    reads_counts_percent = reads_counts_sum.div(total_reads)
    reads_counts_percent = reads_counts_percent * 100

    # Save output tables
    reads_counts_sum.to_csv(
        DATASET_NAME + ".reads.counts.by.organism.table",
        sep="\t",
    )
    reads_counts_percent.to_csv(
        DATASET_NAME + ".reads.percent.by.organism.table",
        sep="\t",
    )
    coverage_counts_sum.to_csv(
        DATASET_NAME + ".coverage.counts.by.organism.table",
        sep="\t",
    )
