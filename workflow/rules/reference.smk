rule reference_recompress_human:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=features["reference"]["human"],
    output:
        fa_gz=REFERENCE / "human.fa.gz",
    log:
        REFERENCE / "human.log",
    conda:
        "../envs/reference.yml"
    threads: 8
    shell:
        """
        (gzip \
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


rule reference_recompress_chicken:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        fa_gz=features["reference"]["chicken"],
    output:
        fa_gz=REFERENCE / "chicken.fa.gz",
    log:
        REFERENCE / "chicken.log",
    conda:
        "../envs/reference.yml"
    threads: 8
    shell:
        """
        (gzip \
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


rule reference_recompress_mag_catalogue_one:
    """Extract the fasta.gz on config.yaml into genome.fa,gz with bgzip"""
    input:
        # fa_gz=features["mag_catalogues"]["{catalogue}"],
        fa_gz=lambda wildcards: features["mag_catalogues"][wildcards.catalogue],
    output:
        fa_gz=REFERENCE / "mags/{catalogue}.fa.gz",
    log:
        REFERENCE / "mags/{catalogue}.log",
    conda:
        "../envs/reference.yml"
    threads: 8
    shell:
        """
        (gzip \
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
    input:
        [REFERENCE / f"mags/{catalogue}.fa.gz" for catalogue in MAG_CATALOGUES],


rule reference:
    """Re-bgzip the reference genome and known variants"""
    input:
        rules.reference_recompress_human.output,
        rules.reference_recompress_chicken.output,
        rules.reference_recompress_mag_catalogue_all.input,
