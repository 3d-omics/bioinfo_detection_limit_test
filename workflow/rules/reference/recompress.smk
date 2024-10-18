rule reference__recompress__mags__:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=lambda wildcards: features["mag_catalogues"][wildcards.catalogue],
    output:
        fa_gz=REFERENCE / "mags" / "{catalogue}.fa.gz",
    log:
        REFERENCE / "mags" / "{catalogue}.log",
    conda:
        "__environment__.yml"
    cache: "omit-software"
    shell:
        """
        ( gzip \
            --decompress \
            --stdout {input.fa_gz} \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
            --stdout \
            /dev/stdin \
        > {output.fa_gz} \
        ) 2> {log}
        """


rule reference__recompress__mags:
    """Recompress all MAG catalogues"""
    input:
        [REFERENCE / "mags" / f"{catalogue}.fa.gz" for catalogue in MAG_CATALOGUES],


rule reference__recompress:
    """Recompress all reference genomes and MAG catalogues"""
    input:
        rules.reference__recompress__mags.input,
