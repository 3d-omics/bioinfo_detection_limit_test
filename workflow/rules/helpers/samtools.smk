rule _helpers__samtools__flagstats_cram:
    """Compute flagstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        txt="{prefix}.flagstats.txt",
    log:
        "{prefix}.flagstats.log",
    conda:
        "_env.yml"
    resources:
        mem_mb=4 * 1024,
    shell:
        "samtools flagstats {input.cram} > {output.txt} 2> {log}"


rule _helpers__samtools__idxstats_cram:
    """Compute idxstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        tsv="{prefix}.idxstats.tsv",
    log:
        "{prefix}.idxstats.log",
    conda:
        "_env.yml"
    resources:
        mem_mb=4 * 1024,
    shell:
        "samtools idxstats {input.cram} > {output.tsv} 2> {log}"


rule _helpers__samtools__stats_cram_host:
    """Compute stats for a host cram"""
    input:
        cram=BOWTIE2_HOSTS / "{genome}" / "{sample}.{library}.cram",
        crai=BOWTIE2_HOSTS / "{genome}" / "{sample}.{library}.cram.crai",
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        tsv=BOWTIE2_HOSTS / "{genome}" / "{sample}.{library}.stats.tsv",
    log:
        BOWTIE2_HOSTS / "{genome}" / "{sample}.{library}.stats.log",
    conda:
        "_env.yml"
    shell:
        """
        samtools stats \
            --reference {input.reference} \
            {input.cram} \
        > {output.tsv} \
        2> {log}
        """


rule _helpers__samtools__stats__cram_mag:
    """Compute stats for a MAG cram"""
    input:
        cram=BOWTIE2_MAGS / "{mag_catalogue}" / "{sample}.{library}.cram",
        crai=BOWTIE2_MAGS / "{mag_catalogue}" / "{sample}.{library}.cram.crai",
        reference=REFERENCE / "mags" / "{mag_catalogue}.fa.gz",
    output:
        tsv=BOWTIE2_MAGS / "{mag_catalogue}" / "{sample}.{library}.stats.tsv",
    log:
        BOWTIE2_MAGS / "{mag_catalogue}" / "{sample}.{library}.stats.log",
    conda:
        "_env.yml"
    shell:
        """
        samtools stats \
            --reference {input.reference} \
            {input.cram} \
        > {output.tsv} \
        2> {log}
        """
