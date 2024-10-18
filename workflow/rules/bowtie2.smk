include: "bowtie2_functions.smk"


rule quantify__bowtie2__build:
    """Build bowtie2 index for the MAG reference"""
    input:
        ref=REFERENCE / "mags" / "{mag_catalogue}.fa.gz",
    output:
        multiext(
            str(QUANT_INDEX) + "/{mag_catalogue}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        QUANT_INDEX / "{mag_catalogue}.log",
    retries: 5
    cache: "omit-software"
    wrapper:
        "v4.7.2/bio/bowtie2/build"


rule quantify__bowtie2__build__all:
    """Build bowtie2 indexes for all the MAG catalogues"""
    input:
        [
            QUANT_INDEX / f"{mag_catalogue}.{extension}"
            for extension in [
                "1.bt2l",
                "2.bt2l",
                "3.bt2l",
                "4.bt2l",
                "rev.1.bt2l",
                "rev.2.bt2l",
            ]
            for mag_catalogue in MAG_CATALOGUES
        ],


rule quantify__bowtie2__map:
    """Map one library to reference genome using bowtie2"""
    input:
        forward_=PRE_BOWTIE2 / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=PRE_BOWTIE2 / "{sample_id}.{library_id}_2.fq.gz",
        mock=multiext(
            str(QUANT_INDEX) + "/{mag_catalogue}",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
        reference=REFERENCE / "mags" / "{mag_catalogue}.fa.gz",
        fai=REFERENCE / "mags" / "{mag_catalogue}.fa.gz.fai",
        gzi=REFERENCE / "mags" / "{mag_catalogue}.fa.gz.gzi",
    output:
        bam=QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.bam",
    log:
        QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.log",
    conda:
        "../environments/bowtie2.yml"
    params:
        samtools_mem=params["quantify"]["bowtie2"]["samtools_mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
        index=lambda w: QUANT_INDEX / f"{w.mag_catalogue}",
    retries: 5
    shell:
        """
        ( bowtie2 \
            -x {params.index} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
            --rg '{params.rg_extra}' \
            --rg-id '{params.rg_id}' \
            --threads {threads} \
        | samtools sort \
            --output-fmt BAM \
            --reference {input.reference} \
            --threads {threads} \
            -T {output.bam} \
            -m {params.samtools_mem} \
            -o {output.bam} \
        ) 2> {log} 1>&2
        """


rule quantify__bowtie2__map__all:
    """Run bowtie2 over all mag catalogues and samples"""
    input:
        [
            QUANT_BOWTIE2 / mag_catalogue / f"{sample_id}.{library_id}.bam"
            for sample_id, library_id in SAMPLE_LIBRARY
            for mag_catalogue in MAG_CATALOGUES
        ],


rule quantify__bowtie2__all:
    """Run bowtie2 build and mappings of all MAG catalogues vs all samples"""
    input:
        rules.quantify__bowtie2__build__all.input,
        rules.quantify__bowtie2__map__all.input,
