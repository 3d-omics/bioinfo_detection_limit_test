rule quantify__bowtie2__cram_to_fastq__:
    input:
        cram=get_host_clean_cram,
    output:
        forward_=temp(
            QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}_1.fq.gz"
        ),
        reverse_=temp(
            QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}_2.fq.gz"
        ),
    log:
        QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.cram_to_fastq.log",
    conda:
        "__environment__.yml"
    shell:
        """
        rm -rf {output.forward_}.collate

        ( samtools view \
            -f 12 \
            -u \
            --threads {threads} \
            {input} \
        | samtools collate \
            -O \
            -u \
            -T {output.forward_}.collate \
            --threads {threads} \
            - \
        | samtools fastq \
            -1 {output.forward_} \
            -2 {output.reverse_} \
            --threads {threads} \
            -c 0 \
            /dev/stdin \
        ) 2> {log} 1>&2
        """


# bowtie2 does not like pipes and/or bams


rule quantify__bowtie2__:
    """Map one library to reference genome using bowtie2"""
    input:
        forward_=QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}_2.fq.gz",
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
        ( bowtie2 \
            -x {params.index} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
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
