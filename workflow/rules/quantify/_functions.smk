def get_coverm_genome_tsv_files_for_aggregation(wildcards):
    """Get all coverm genome tsv files for aggregation"""
    mag_catalogue = wildcards.mag_catalogue
    method = wildcards.method
    return [
        COVERM / mag_catalogue / "genome" / method / f"{sample_id}.{library_id}.tsv"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]


def get_coverm_contig_tsv_files_for_aggregation(wildcards):
    """Get all coverm genome tsv files for aggregation"""
    mag_catalogue = wildcards.mag_catalogue
    method = wildcards.method
    return [
        COVERM / mag_catalogue / "contig" / method / f"{sample_id}.{library_id}.tsv"
        for sample_id, library_id in SAMPLE_LIBRARY
    ]
