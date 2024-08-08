def get_reads_reports_for_library_reports(wildcards):
    """Compose the paths for the reads reports"""
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    return [READS / f"{sample_id}.{library_id}_{end}_fastqc.zip" for end in ["1", "2"]]


def get_fastp_reports_for_library_reports(wildcards):
    """Compose the paths for the fastp reports"""
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    return [
        FASTP / f"{sample_id}.{library_id}_fastp.json",
        FASTP / f"{sample_id}.{library_id}_1_fastqc.zip",
        FASTP / f"{sample_id}.{library_id}_2_fastqc.zip",
    ]


def get_nonpareil_reports_for_library_reports(wildcards):
    """Compose the paths to the nonpareil json reports"""
    return NONPAREIL / f"{wildcards.sample_id}.{wildcards.library_id}.json"


def get_bowtie2_host_for_library_reports(wildcards):
    """Compose the paths for the bowtie2_hosts reports"""
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    return [
        PRE_BOWTIE2 / host_name / f"{sample_id}.{library_id}.{report}"
        for host_name in HOST_NAMES
        for report in BAM_REPORTS
    ]


def get_kraken2_for_library_reports(wildcards):
    """Compose the paths for the kraken2 reports"""
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    return [
        KRAKEN2 / kraken2_db / f"{sample_id}.{library_id}.report"
        for kraken2_db in KRAKEN2_DBS
    ]


def get_bowtie2_mags_for_library_reports(wildcards):
    """Compose the paths for the bowtie2_mags reports"""
    sample_id = wildcards.sample_id
    library_id = wildcards.library_id
    return [
        QUANT_BOWTIE2 / mag_catalogue / f"{sample_id}.{library_id}.{report}"
        for mag_catalogue in MAG_CATALOGUES
        for report in ["stats.tsv", "flagstats.txt"]
    ]
