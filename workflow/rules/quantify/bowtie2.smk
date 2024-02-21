rule quantify__bowtie2__:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        cram=get_host_clean_cram,
        mock=multiext(
            str(QUANT_INDEX) + "/{mag_catalogue}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
        reference=REFERENCE / "mags" / "{mag_catalogue}.fa.gz",
    output:
        cram=QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.cram",
    log:
        QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    threads: 24
    resources:
        mem_mb=double_ram(params["quantify"]["bowtie2"]["mem_gb"]),
        runtime=24 * 60,
    params:
        samtools_mem=params["quantify"]["bowtie2"]["samtools_mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
        index=lambda w: QUANT_INDEX / f"{w.mag_catalogue}",
    group:
        "sample"
    retries: 5
    shell:
        """
        ( samtools view \
            -u \
            -f 12 \
            {input.cram} \
        | samtools sort \
            -u \
            -n \
            --threads {threads} \
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
