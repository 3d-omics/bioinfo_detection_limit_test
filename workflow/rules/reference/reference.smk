rule reference_recompress_host_one:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=lambda wildcards: features["references"][wildcards.genome],
    output:
        fa_gz=REFERENCE / "{genome}.fa.gz",
    log:
        REFERENCE / "{genome}.log",
    conda:
        "reference.yml"
    threads: 8
    shell:
        """
        ( gzip \
            --decompres \
            --stdout {input.fa_gz} \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
            --stdout \
            /dev/stdin \
        > {output.fa_gz} \
        ) 2> {log}
        """


rule reference_recompress_host_all:
    """Recompress all host genomes"""
    input:
        [REFERENCE / f"{genome}.fa.gz" for genome in HOST_NAMES],


rule reference_recompress_mag_catalogue_one:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=lambda wildcards: features["mag_catalogues"][wildcards.catalogue],
    output:
        fa_gz=REFERENCE / "mags" / "{catalogue}.fa.gz",
    log:
        REFERENCE / "mags" / "{catalogue}.log",
    conda:
        "reference.yml"
    threads: 8
    shell:
        """
        ( gzip \
            --decompres \
            --stdout {input.fa_gz} \
        | bgzip \
            --compress-level 9 \
            --threads {threads} \
            --stdout \
            /dev/stdin \
        > {output.fa_gz} \
        ) 2> {log}
        """


rule reference_recompress_mag_catalogue_all:
    """Recompress all MAG catalogues"""
    input:
        [REFERENCE / "mags" / f"{catalogue}.fa.gz" for catalogue in MAG_CATALOGUES],


rule reference:
    """Re-bgzip the reference genome and known variants"""
    input:
        rules.reference_recompress_host_one.input,
        rules.reference_recompress_mag_catalogue_all.input,
