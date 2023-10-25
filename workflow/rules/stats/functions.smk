def get_input_forward_for_stats(wildcards):
    """Get the forward input file for stats"""
    last_genome = HOST_NAMES[-1]
    if len(HOST_NAMES) == 0:
        return FASTP / f"{wildcards.sample}.{wildcards.library}_1.fq.gz"
    return (
        BOWTIE2_HOSTS
        / f"non{last_genome}/{wildcards.sample}.{wildcards.library}_1.fq.gz"
    )


def get_input_reverse_for_stats(wildcards):
    """Get the forward input file for stats"""
    last_genome = HOST_NAMES[-1]
    if len(HOST_NAMES) == 0:
        return FASTP / f"{wildcards.sample}.{wildcards.library}_2.fq.gz"
    return (
        BOWTIE2_HOSTS
        / f"non{last_genome}/{wildcards.sample}.{wildcards.library}_2.fq.gz"
    )


def compose_prefix_for_nonpareil(wildcards):
    """Compose prefix for nonpareil output files"""
    return STATS_NONPAREIL / f"{wildcards.sample}.{wildcards.library}"


def get_coverm_genome_tsv_files_for_aggregation(wildcards):
    """Get all coverm genome tsv files for aggregation"""
    mag_catalogue = wildcards.mag_catalogue
    method = wildcards.method
    return [
        STATS_COVERM / f"{mag_catalogue}/genome/{method}/{sample}.{library}.tsv"
        for sample, library in SAMPLE_LIB
    ]


def get_coverm_contig_tsv_files_for_aggregation(wildcards):
    """Get all coverm genome tsv files for aggregation"""
    mag_catalogue = wildcards.mag_catalogue
    method = wildcards.method
    return [
        STATS_COVERM / f"{mag_catalogue}/contig/{method}/{sample}.{library}.tsv"
        for sample, library in SAMPLE_LIB
    ]


def get_kraken2_database(wildcards):
    return features["kraken2_dbs"][wildcards.kraken2_db]


def compose_out_folder_for_pre_kraken2_assign_all(wildcards):
    return KRAKEN2 / f"{wildcards.kraken2_db}"
