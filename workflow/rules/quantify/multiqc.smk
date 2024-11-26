rule quantify__multiqc:
    """Collect all reports for the bowtie2 step when mapping to a mag catalogue"""
    input:
        bowtie2=[
            QUANT_BOWTIE2 / f"{mag_catalogue}.{sample_id}.{library_id}.stats.tsv"
            for mag_catalogue in MAG_CATALOGUES
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        RESULTS / "quantify.html",
        RESULTS / "quantify.zip",
    log:
        RESULTS / "quantify.log",
    params:
        extra="--title quantify --dirs --dirs-depth 1 --fullnames --force",
    resources:
        mem_mb=double_ram(8 * 1024),
        runtime=6 * 60,
    wrapper:
        "v5.1.0/bio/multiqc"


rule quantify__multiqc__all:
    """Collect all per step reports for the pipeline"""
    input:
        rules.quantify__multiqc.output,
