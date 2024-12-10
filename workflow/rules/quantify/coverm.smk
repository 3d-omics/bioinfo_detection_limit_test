# coverm genome ----
use rule coverm__genome as quantify__coverm__genome with:
    input:
        QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.bam",
    output:
        temp(
            QUANT_COVERM
            / "{mag_catalogue}"
            / "genome.{method}.{sample_id}.{library_id}.tsv.gz"
        ),
    log:
        QUANT_COVERM
        / "{mag_catalogue}"
        / "genome.{method}.{sample_id}.{library_id}.log",
    params:
        method=lambda w: w.method,
        extra=params["quantify"]["coverm"]["genome"]["extra"],
        separator=params["quantify"]["coverm"]["genome"]["separator"],


rule quantify__coverm__genome__join:
    input:
        lambda w: [
            QUANT_COVERM
            / f"{w.mag_catalogue}"
            / f"genome.{w.method}.{sample_id}.{library_id}.tsv.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        QUANT_COVERM / "genome.{method, \w+}.{mag_catalogue}.tsv.gz",
    log:
        QUANT_COVERM / "genome.{method, \w+}.{mag_catalogue}.log",
    params:
        subcommand="join",
    wrapper:
        "v5.2.1/utils/csvtk"


rule quantify__coverm__genome__all:
    """Run coverm genome and all methods"""
    input:
        [
            QUANT_COVERM / f"genome.{method}.{mag_catalogue}.tsv.gz"
            for method in params["quantify"]["coverm"]["genome"]["methods"]
            for mag_catalogue in MAG_CATALOGUES
        ],


rule quantify__coverm__all:
    """Run both coverm overall and contig"""
    input:
        rules.quantify__coverm__genome__all.input,
