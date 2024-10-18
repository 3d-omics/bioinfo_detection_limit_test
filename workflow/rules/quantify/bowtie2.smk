rule quantify__bowtie2__:
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
        "__environment__.yml"
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


rule quantify__bowtie2:
    """Run bowtie2 over all mag catalogues and samples"""
    input:
        [
            QUANT_BOWTIE2 / mag_catalogue / f"{sample_id}.{library_id}.bam"
            for sample_id, library_id in SAMPLE_LIBRARY
            for mag_catalogue in MAG_CATALOGUES
        ],
