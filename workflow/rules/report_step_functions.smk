def get_report_step_kraken2_reports(wildcards):
    """Get all reports for the kraken2 step"""
    kraken2_db = wildcards.kraken2_db
    return [
        KRAKEN2 / f"{kraken2_db}/{sample}.{library}.report"
        for sample, library in SAMPLE_LIB
    ]
