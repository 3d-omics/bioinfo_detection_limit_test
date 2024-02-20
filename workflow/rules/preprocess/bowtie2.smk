rule _preprocess__bowtie2__build:
    """Build bowtie2 index for the human reference

    Let the script decide to use a small or a large index based on the size of
    the reference genome.
    """
    input:
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        mock=touch(PRE_BOWTIE2 / "{genome}_index"),
    log:
        PRE_BOWTIE2 / "{genome}_index.log",
    conda:
        "__environment__.yml"
    threads: 24
    resources:
        mem_mb=double_ram(params["preprocess"]["bowtie2"]["mem_gb"]),
        runtime=24 * 60,
    retries: 5
    cache: True
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.reference} \
            {output.mock} \
        2> {log} 1>&2
        """


rule _preprocess__bowtie2__map:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        cram=get_input_cram_for_host_mapping,
        mock=PRE_BOWTIE2 / "{genome}_index",
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        cram=PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.cram",
    log:
        PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.log",
    params:
        samtools_mem=params["preprocess"]["bowtie2"]["samtools_mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    threads: 24
    conda:
        "__environment__.yml"
    resources:
        mem_mb=double_ram(params["preprocess"]["bowtie2"]["mem_gb"]),
        runtime=24 * 60,
    group:
        "sample"
    retries: 5
    shell:
        """
        ( samtools view \
            -f 12 \
            -u \
            {input.cram} \
        | samtools sort \
            -u \
            -n \
            --threads {threads} \
        | bowtie2 \
            -x {input.mock} \
            -b /dev/stdin \
            --align-paired-reads \
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
        | samtools sort \
            --output-fmt CRAM \
            --reference {input.reference} \
            --threads {threads} \
            -T {output.cram} \
            -m {params.samtools_mem} \
            -o {output.cram}
        ) 2> {log} 1>&2
        """


rule preprocess__bowtie2:
    """Run bowtie2 on all libraries and generate reports"""
    input:
        [
            PRE_BOWTIE2 / "non_host" / f"{sample_id}.{library_id}_{end}.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in [1, 2]
        ],
