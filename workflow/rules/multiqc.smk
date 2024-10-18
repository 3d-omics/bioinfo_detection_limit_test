rule quantify__multiqc:
    """Collect all reports for the bowtie2 step when mapping to a mag catalogue"""
    input:
        reports=[
            QUANT_BOWTIE2 / mag_catalogue / f"{sample_id}.{library_id}.{report}"
            for mag_catalogue in MAG_CATALOGUES
            for sample_id, library_id in SAMPLE_LIBRARY
            for report in ["stats.tsv", "flagstats.txt"]
        ],
    output:
        RESULTS / "quantify.html",
    log:
        RESULTS / "quantify.log",
    conda:
        "../environments/multiqc.yml"
    params:
        dir=RESULTS,
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


rule quantify__multiqc__all:
    """Collect all per step reports for the pipeline"""
    input:
        rules.quantify__multiqc.output,
