rule _reference__recompress_host:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=lambda wildcards: features["references"][wildcards.genome],
    output:
        fa_gz=REFERENCE / "{genome}.fa.gz",
    log:
        REFERENCE / "{genome}.log",
    conda:
        "_env.yml"
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


rule reference__recompress_hosts:
    """Recompress all host genomes"""
    input:
        [REFERENCE / f"{genome}.fa.gz" for genome in HOST_NAMES],


rule _reference__recompress_mags:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=lambda wildcards: features["mag_catalogues"][wildcards.catalogue],
    output:
        fa_gz=REFERENCE / "mags" / "{catalogue}.fa.gz",
    log:
        REFERENCE / "mags" / "{catalogue}.log",
    conda:
        "_env.yml"
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


rule reference__recompress_mags:
    """Recompress all MAG catalogues"""
    input:
        [REFERENCE / "mags" / f"{catalogue}.fa.gz" for catalogue in MAG_CATALOGUES],


rule reference:
    """Re-bgzip the reference genome and known variants"""
    input:
        rules.reference__recompress_hosts.input,
        rules.reference__recompress_mags.input,
