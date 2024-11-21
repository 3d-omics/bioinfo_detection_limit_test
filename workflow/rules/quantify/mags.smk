rule quantify__mags:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        lambda w: features["mag_catalogues"][w.catalogue],
    output:
        QUANT_MAGS / "{catalogue}.fa.gz",
    log:
        QUANT_MAGS / "{catalogue}.log",
    conda:
        "../../environments/mags.yml"
    cache: "omit-software"
    shell:
        """
        ( gzip \
            --decompress \
            --stdout {input} \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
            --stdout \
            /dev/stdin \
        > {output} \
        ) 2> {log}
        """


rule quantify__mags__all:
    """Recompress all MAG catalogues"""
    input:
        [QUANT_MAGS / f"{catalogue}.fa.gz" for catalogue in MAG_CATALOGUES],
