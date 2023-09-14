# rule crai:
#     """Generate a cram index"""
#     input:
#         cram="{prefix}.cram",
#     output:
#         crai="{prefix}.cram.crai",
#     log:
#         "{prefix}.cram.crai.log",
#     conda:
#         "../envs/samtools.yml"
#     shell:
#         "samtools index {input} 2> {log} 1>&2"


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
        cram=BOWTIE2_HOSTS / "{genome}/{sample}.{library}_{library_type}.cram",
        crai=BOWTIE2_HOSTS / "{genome}/{sample}.{library}_{library_type}.cram.crai",
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        tsv=BOWTIE2_HOSTS / "{genome}/{sample}.{library}_{library_type}.stats.tsv",
    log:
        BOWTIE2_HOSTS / "{genome}/{sample}.{library}_{library_type}.stats.log",
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
        cram=BOWTIE2_MAGS / "{mag_catalogue}/{sample}.{library}_{library_type}.cram",
        crai=BOWTIE2_MAGS
        / "{mag_catalogue}/{sample}.{library}_{library_type}.cram.crai",
        reference=REFERENCE / "mags/{mag_catalogue}.fa.gz",
    output:
        tsv=BOWTIE2_MAGS / "{mag_catalogue}/{sample}.{library}_{library_type}.stats.tsv",
    log:
        BOWTIE2_MAGS / "{mag_catalogue}/{sample}.{library}_{library_type}.stats.log",
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
