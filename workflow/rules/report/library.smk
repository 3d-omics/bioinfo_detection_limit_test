rule report__library__:
    """Make a MultiQC report for a single library"""
    input:
        get_reads_reports_for_library_reports,
        get_fastp_reports_for_library_reports,
        get_nonpareil_reports_for_library_reports,
        get_bowtie2_host_for_library_reports,
        get_kraken2_for_library_reports,
        get_bowtie2_mags_for_library_reports,
    output:
        REPORT_LIBRARY / "{sample_id}.{library_id}.html",
    log:
        REPORT_LIBRARY / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    params:
        library=lambda w: f"{w.sample_id}.{w.library_id}",
        out_dir=REPORT_LIBRARY,
    retries: 5
    shell:
        """
        multiqc \
            --title {params.library} \
            --force \
            --filename {params.library} \
            --outdir {params.out_dir} \
            --dirs \
            --dirs-depth 1 \
            {input} \
        2> {log} 1>&2
        """


rule report__library:
    """Make a MultiQC report for every library"""
    input:
        [
            REPORT_LIBRARY / f"{sample_id}.{library_id}.html"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
