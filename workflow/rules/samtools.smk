rule samtools_flagstats_cram:
    """Compute flagstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        txt="{prefix}.flagstats.txt",
    log:
        "{prefix}.flagstats.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools flagstats {input.cram} > {output.txt} 2> {log}"


rule samtools_idxstats_cram:
    """Compute idxstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        tsv="{prefix}.idxstats.tsv",
    log:
        "{prefix}.idxstats.log",
    conda:
        "../envs/samtools.yml"
    shell:
        "samtools idxstats {input.cram} > {output.tsv} 2> {log}"


rule samtools_stats_cram_host:
    """Compute stats for a cram"""
    input:
        cram=BOWTIE2_HOSTS / "{genome}/{sample}.{library}.cram",
        crai=BOWTIE2_HOSTS / "{genome}/{sample}.{library}.cram.crai",
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        tsv=BOWTIE2_HOSTS / "{genome}/{sample}.{library}.stats.tsv",
    log:
        BOWTIE2_HOSTS / "{genome}/{sample}.{library}.stats.log",
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools stats \
            --reference {input.reference} \
            {input.cram} \
        > {output.tsv} \
        2> {log}
        """


rule samtools_stats_cram_mag_catalogue:
    input:
        cram=BOWTIE2_MAGS / "{mag_catalogue}/{sample}.{library}.cram",
        crai=BOWTIE2_MAGS / "{mag_catalogue}/{sample}.{library}.cram.crai",
        reference=REFERENCE / "mags/{mag_catalogue}.fa.gz",
    output:
        tsv=BOWTIE2_MAGS / "{mag_catalogue}/{sample}.{library}.stats.tsv",
    log:
        BOWTIE2_MAGS / "{mag_catalogue}/{sample}.{library}.stats.log",
    conda:
        "../envs/samtools.yml"
    shell:
        """
        samtools stats \
            --reference {input.reference} \
            {input.cram} \
        > {output.tsv} \
        2> {log}
        """
