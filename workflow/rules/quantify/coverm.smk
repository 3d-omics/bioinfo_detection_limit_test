# coverm genome ----
use rule coverm__genome as quantify__coverm__genome with:
    input:
        QUANT_BOWTIE2 / "{catalogue}.{sample_id}.{library_id}.bam",
    output:
        temp(
            QUANT_COVERM
            / "files"
            / "genome.{method}.{catalogue}.{sample_id}.{library_id}.tsv.gz"
        ),
    log:
        QUANT_COVERM
        / "files"
        / "genome.{method}.{catalogue}.{sample_id}.{library_id}.log",
    params:
        method=lambda w: w.method,
        extra=params["quantify"]["coverm"]["genome"]["extra"],
        separator=params["quantify"]["coverm"]["genome"]["separator"],


use rule csvkit__csvjoin as quantify__coverm__genome__csvjoin with:
    input:
        lambda w: [
            QUANT_COVERM
            / "files"
            / f"genome.{w.method}.{w.catalogue}.{sample_id}.{library_id}.tsv.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        QUANT_COVERM / "coverm.genome.{method}.{catalogue}.tsv.gz",
    log:
        QUANT_COVERM / "coverm.genome.{method}.{catalogue}.log",
    conda:
        "../../environments/csvkit.yml"


rule quantify__coverm__genome__all:
    """Run coverm genome and all methods"""
    input:
        [
            QUANT_COVERM / f"coverm.genome.{method}.{catalogue}.tsv.gz"
            for method in params["quantify"]["coverm"]["genome"]["methods"]
            for catalogue in MAG_CATALOGUES
        ],


# coverm contig ----
use rule coverm__contig as quantify__coverm__contig with:
    input:
        QUANT_BOWTIE2 / "{catalogue}.{sample_id}.{library_id}.bam",
    output:
        temp(
            QUANT_COVERM
            / "files"
            / "contig.{method}.{catalogue}.{sample_id}.{library_id}.tsv.gz"
        ),
    log:
        QUANT_COVERM
        / "files"
        / "contig.{method}.{catalogue}.{sample_id}.{library_id}.log",


use rule csvkit__csvjoin as quantify__coverm__contig__csvjoin with:
    input:
        lambda w: [
            QUANT_COVERM
            / "files"
            / f"contig.{w.method}.{w.mag_catalogue}.{sample_id}.{library_id}.tsv.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    output:
        QUANT_COVERM / "coverm.contig.{method}.{mag_catalogue}.tsv.gz",
    log:
        QUANT_COVERM / "coverm.contig.{method}.{mag_catalogue}.log",
    conda:
        "../../environments/csvkit.yml"


rule quantify__coverm__contig__all:
    """Run all rules to run coverm contig over all MAG catalogues"""
    input:
        [
            QUANT_COVERM / f"coverm.contig.{method}.{mag_catalogue}.tsv.gz"
            for mag_catalogue in MAG_CATALOGUES
            for method in params["quantify"]["coverm"]["contig"]["methods"]
        ],


rule quantify__coverm__all:
    """Run both coverm overall and contig"""
    input:
        rules.quantify__coverm__genome__all.input,
        rules.quantify__coverm__contig__all.input,
