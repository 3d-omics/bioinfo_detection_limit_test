def get_bowtie2_host_for_library_reports(wildcards):
    """Compose the paths for the bowtie2_hosts reports"""
    sample = wildcards.sample
    library = wildcards.library
    library_type = "pe" if [sample, library] in SAMPLE_LIB_PE else "se"
    return [
        BOWTIE2_HOSTS / f"{host_name}/{sample}.{library}_{library_type}.{report}"
        for host_name in HOST_NAMES
        for report in BAM_REPORTS
    ]


def get_bowtie2_mags_for_library_reports(wildcards):
    """Compose the paths for the bowtie2_mags reports"""
    sample = wildcards.sample
    library = wildcards.library
    library_type = "pe" if [sample, library] in SAMPLE_LIB_PE else "se"
    return [
        BOWTIE2_MAGS / f"{mag_catalogue}/{sample}.{library}_{library_type}.{report}"
        for mag_catalogue in MAG_CATALOGUES
        for report in BAM_REPORTS
    ]


def get_kraken2_for_library_reports(wildcards):
    """Compose the paths for the kraken2 reports"""
    sample = wildcards.sample
    library = wildcards.library
    library_type = "pe" if [sample, library] in SAMPLE_LIB_PE else "se"
    return [
        f"{KRAKEN2}/{kraken2_db}/{sample}.{library}_{library_type}.report"
        for kraken2_db in KRAKEN2_DBS
    ]


def get_report_step_kraken2_reports(wildcards):
    """Get all reports for the kraken2 step"""
    kraken2_db = wildcards.kraken2_db
    return [
        KRAKEN2 / f"{kraken2_db}/{sample}.{library}_pe.report"
        for sample, library in SAMPLE_LIB_PE
    ] + [
        KRAKEN2 / f"{kraken2_db}/{sample}.{library}_se.report"
        for sample, library in SAMPLE_LIB_SE
    ]
