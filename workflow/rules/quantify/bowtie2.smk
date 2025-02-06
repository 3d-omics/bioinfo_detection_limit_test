use rule bowtie2__build from mg_assembly as quantify__bowtie2__build with:
    input:
        ref=QUANT_MAGS / "{mag_catalogue}.fa.gz",
    output:
        multiext(
            str(QUANT_BUILD / "{mag_catalogue}."),
            "1.bt2",
            "2.bt2",
            "3.bt2",
            "4.bt2",
            "rev.1.bt2",
            "rev.2.bt2",
        ),
    log:
        QUANT_BUILD / "{mag_catalogue}.log",
    cache: "omit-software"
    retries: 5


rule quantify__bowtie2__build__all:
    input:
        [
            QUANT_BUILD / f"{mag_catalogue}.{extension}"
            for mag_catalogue in MAG_CATALOGUES
            for extension in [
                "1.bt2",
                "2.bt2",
                "3.bt2",
                "4.bt2",
                "rev.1.bt2",
                "rev.2.bt2",
            ]
        ],


use rule bowtie2__map from mg_assembly as quantify__bowtie2__map with:
    input:
        forward_=PRE_CLEAN / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_CLEAN / "{sample_id}.{library_id}_2.fq.gz",
        mock=multiext(
            str(QUANT_BUILD / "{mag_catalogue}"),
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    output:
        QUANT_BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.bam",
    log:
        QUANT_BOWTIE2 / "{mag_catalogue}.{sample_id}.{library_id}.log",
    params:
        index=lambda w: QUANT_BUILD / f"{w.mag_catalogue}",
        bowtie2_extra=params["quantify"]["bowtie2"]["bowtie2_extra"],
        samtools_extra=params["quantify"]["bowtie2"]["samtools_extra"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    retries: 5


rule quantify__bowtie2__map__all:
    """Run bowtie2 over all mag catalogues and samples"""
    input:
        [
            QUANT_BOWTIE2 / f"{mag_catalogue}.{sample_id}.{library_id}.bam"
            for sample_id, library_id in SAMPLE_LIBRARY
            for mag_catalogue in MAG_CATALOGUES
        ],


rule quantify__bowtie2__all:
    """Run bowtie2 build and mappings of all MAG catalogues vs all samples"""
    input:
        rules.quantify__bowtie2__build__all.input,
        rules.quantify__bowtie2__map__all.input,
