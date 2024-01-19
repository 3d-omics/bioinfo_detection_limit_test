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
