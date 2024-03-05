rule quantify__bowtie2__:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        cram=get_host_clean_cram,
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
        cram=QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.cram",
    log:
        QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    params:
        samtools_mem=params["quantify"]["bowtie2"]["samtools_mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
        index=lambda w: QUANT_INDEX / f"{w.mag_catalogue}",
    retries: 5
    shell:
        """
        ( samtools view \
            -u \
            -f 12 \
            {input.cram} \
        | samtools collate \
            -u \
            -O \
            --threads {threads} \
            - \
        | bowtie2 \
            -x {params.index} \
            -b /dev/stdin \
            --align-paired-reads \
            --rg '{params.rg_extra}' \
            --rg-id '{params.rg_id}' \
            --threads {threads} \
        | samtools sort \
            --output-fmt CRAM \
            --reference {input.reference} \
            --threads {threads} \
            -T {output.cram} \
            -m {params.samtools_mem} \
            -o {output.cram} \
        ) 2> {log} 1>&2
        """


rule quantify__bowtie2:
    """Run bowtie2 over all mag catalogues and samples"""
    input:
        [
            QUANT_BOWTIE2 / mag_catalogue / f"{sample_id}.{library_id}.cram"
            for sample_id, library_id in SAMPLE_LIBRARY
            for mag_catalogue in MAG_CATALOGUES
        ],
