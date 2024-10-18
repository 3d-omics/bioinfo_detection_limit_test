def get_coverm_file(wildcards, coverm_type):
    assert coverm_type in ["genome", "contig"]
    mag_catalogue = wildcards.mag_catalogue
    method = wildcards.method
    return [
        COVERM / mag_catalogue / coverm_type / method / f"{sample_id}.{library_id}.tsv"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]


def get_coverm_genome_tsv_files_for_aggregation(wildcards):
    """Get all coverm genome tsv files for aggregation"""
    return get_coverm_file(wildcards, "genome")


def get_coverm_contig_tsv_files_for_aggregation(wildcards):
    """Get all coverm genome tsv files for aggregation"""
    return get_coverm_file(wildcards, "contig")


def compose_rg_id(wildcards):
    """Compose the read group ID for bowtie2"""
    return f"{wildcards.sample_id}_{wildcards.library_id}"


def compose_rg_extra(wildcards):
    """Compose the read group extra information for bowtie2"""
    return f"LB:truseq_{wildcards.library_id}\\tPL:Illumina\\tSM:{wildcards.sample_id}"
