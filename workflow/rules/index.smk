rule quantify__index__:
    """Build bowtie2 index for the MAG reference

    NOTE: the index is built big on purpose. Bowtie2 creates different indexes depending on the size, bt2 for small and bt2l for large, and it is a pain in the ass to separate cases.
    """
    input:
        reference=REFERENCE / "mags" / "{mag_catalogue}.fa.gz",
    output:
        multiext(
            str(QUANT_INDEX) + "/{mag_catalogue}",
            ".1.bt2l",
            ".2.bt2l",
            ".3.bt2l",
            ".4.bt2l",
            ".rev.1.bt2l",
            ".rev.2.bt2l",
        ),
    log:
        QUANT_INDEX / "{mag_catalogue}.log",
    conda:
        "__environment__.yml"
    params:
        prefix=lambda w: str(QUANT_INDEX / f"{w.mag_catalogue}"),
    retries: 5
    cache: "omit-software"
    shell:
        """
        bowtie2-build \
            --large-index \
            --threads {threads} \
            {input.reference} \
            {params.prefix} \
        2> {log} 1>&2
        """


rule quantify__index:
    """Build bowtie2 indexes for all the MAG catalogues"""
    input:
        [
            QUANT_INDEX / f"{mag_catalogue}.{end}"
            for end in [
                "1.bt2l",
                "2.bt2l",
                "3.bt2l",
                "4.bt2l",
                "rev.1.bt2l",
                "rev.2.bt2l",
            ]
            for mag_catalogue in MAG_CATALOGUES
        ],
