rule preprocess__index__:
    """Build bowtie2 index for a reference

    Let the script decide to use a small or a large index based on the size of
    the reference genome.
    """
    input:
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        multiext(
            str(PRE_INDEX) + "/{genome}",
            ".1.bt2",
            ".2.bt2",
            ".3.bt2",
            ".4.bt2",
            ".rev.1.bt2",
            ".rev.2.bt2",
        ),
    log:
        PRE_INDEX / "{genome}.log",
    conda:
        "__environment__.yml"
    params:
        prefix=lambda w: PRE_INDEX / w.genome,
    threads: 24
    resources:
        mem_mb=double_ram(params["preprocess"]["bowtie2"]["mem_gb"]),
        runtime=24 * 60,
    # retries: 5
    cache: True
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.reference} \
            {params.prefix} \
        2> {log} 1>&2
        """


rule preprocess__index:
    input:
        [
            PRE_INDEX / f"{genome}.{end}"
            for end in ["1.bt2", "2.bt2", "3.bt2", "4.bt2", "rev.1.bt2", "rev.2.bt2"]
            for genome in HOST_NAMES
        ],
