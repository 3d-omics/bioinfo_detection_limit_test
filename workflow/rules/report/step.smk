rule report__step__quantify:
    """Collect all reports for the bowtie2 step when mapping to a mag catalogue"""
    input:
        reports=[
            QUANT_BOWTIE2 / mag_catalogue / f"{sample_id}.{library_id}.{report}"
            for mag_catalogue in MAG_CATALOGUES
            for sample_id, library_id in SAMPLE_LIBRARY
            for report in ["stats.tsv", "flagstats.txt"]
        ],
    output:
        REPORT_STEP / "quantify.html",
    log:
        REPORT_STEP / "quantify.log",
    conda:
        "__environment__.yml"
    params:
        dir=REPORT_STEP,
    retries: 5
    shell:
        """
        multiqc \
            --title quantify \
            --force \
            --filename quantify \
            --outdir {params.dir} \
            --dirs \
            --dirs-depth 1 \
            {input.reports} \
        2> {log} 1>&2
        """


rule report__step:
    """Collect all per step reports for the pipeline"""
    input:
        rules.report__step__quantify.output,
