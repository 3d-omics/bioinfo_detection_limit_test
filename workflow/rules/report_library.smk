rule report_library_one:
    """Make a MultiQC report for a single library"""
    input:
        READS / "{sample}.{library}_1_fastqc.zip",
        READS / "{sample}.{library}_2_fastqc.zip",
        FASTP / "{sample}.{library}_fastp.json",
        FASTP / "{sample}.{library}_1_fastqc.zip",
        FASTP / "{sample}.{library}_2_fastqc.zip",
        get_kraken2_for_library_reports,
        get_bowtie2_host_for_library_reports,
        get_bowtie2_mags_for_library_reports,
    output:
        REPORT_LIBRARY / "{sample}.{library}_{library_type}.html",
    log:
        REPORT_LIBRARY / "{sample}.{library}_{library_type}.log",
    conda:
        "../envs/report.yml"
    params:
        library="{sample}.{library}_{library_type}",
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
            REPORT_LIBRARY / f"{sample}.{library}_pe.html"
            for sample, library in SAMPLE_LIB_PE
        ]
        + [
            REPORT_LIBRARY / f"{sample}.{library}_se.html"
            for sample, library in SAMPLE_LIB_SE
        ],


rule report_library:
    """Make all MultiQC reports per library"""
    input:
        rules.report_library_all.input,


localrules:
    report_library_one,
