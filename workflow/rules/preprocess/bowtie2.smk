rule preprocess__bowtie2__:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        cram=get_input_cram_for_host_mapping,
        mock=multiext(
            str(PRE_INDEX) + "/{genome}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
        reference=REFERENCE / "{genome}.fa.gz",
        fai=REFERENCE / "{genome}.fa.gz.fai",
        gzi=REFERENCE / "{genome}.fa.gz.gzi",
    output:
        cram=PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.cram",
    log:
        PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.log",
    params:
        index=lambda w: PRE_INDEX / f"{w.genome}",
        samtools_mem=params["preprocess"]["bowtie2"]["samtools_mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    conda:
        "__environment__.yml"
    shell:
        """
        ( samtools collate \
            -f \
            -O \
            -T {output.cram}.collate \
            {input.cram} \
        | samtools view \
            -f 12 \
            -1 \
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
        # """
        # ( samtools view \
        #     -u \
        #     -f 12 \
        #     {input.cram} \
        # | samtools sort \
        #     -n \
        #     -u \
        #     -T {output.cram}.sort_by_name \
        #     --threads {threads} \
        #     -m {params.samtools_mem} \
        # | bowtie2 \
        #     -x {params.index} \
        #     -b /dev/stdin \
        #     --align-paired-reads \
        #     --rg '{params.rg_extra}' \
        #     --rg-id '{params.rg_id}' \
        #     --threads {threads} \
        # | samtools sort \
        #     --output-fmt CRAM \
        #     --reference {input.reference} \
        #     --threads {threads} \
        #     -T {output.cram} \
        #     -m {params.samtools_mem} \
        #     -o {output.cram} \
        # ) 2> {log} 1>&2
        # """


rule preprocess__bowtie2:
    """Run bowtie2 on all libraries and generate reports"""
    input:
        [
            PRE_BOWTIE2 / LAST_HOST / f"{sample_id}.{library_id}.cram"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
