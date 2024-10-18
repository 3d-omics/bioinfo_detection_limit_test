rule quantify__mags:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=lambda wildcards: features["mag_catalogues"][wildcards.catalogue],
    output:
        fa_gz=MAGS / "{catalogue}.fa.gz",
    log:
        MAGS / "{catalogue}.log",
    conda:
        "../environments/mags.yml"
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


rule quantify__mags__all:
    """Recompress all MAG catalogues"""
    input:
        [MAGS / f"{catalogue}.fa.gz" for catalogue in MAG_CATALOGUES],
