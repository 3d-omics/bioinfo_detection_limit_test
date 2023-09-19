rule report_library_one:
    """Make a MultiQC report for a single library"""
    input:
        get_reads_reports_for_library_reports,
        get_fastp_reports_for_library_reports,
        get_kraken2_for_library_reports,
        get_bowtie2_host_for_library_reports,
        get_bowtie2_mags_for_library_reports,
    output:
        REPORT_LIBRARY / "{sample}.{library}.html",
    log:
        REPORT_LIBRARY / "{sample}.{library}.log",
    conda:
        "../envs/report.yml"
    params:
        library="{sample}.{library}",
        out_dir=REPORT_LIBRARY,
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


rule report_library_all:
    """Make a MultiQC report for every library"""
    input:
        [
            REPORT_LIBRARY / f"{sample}.{library}.html"
            for sample, library in SAMPLE_LIB_PE
        ],


rule report_library:
    """Make all MultiQC reports per library"""
    input:
        rules.report_library_all.input,


localrules:
    report_library_one,
