rule _helpers__samtools__stats_cram_host:
    """Compute stats for a host cram"""
    input:
        cram=PRE_BOWTIE2 / "{genome}" / "{sample}.{library}.cram",
        crai=PRE_BOWTIE2 / "{genome}" / "{sample}.{library}.cram.crai",
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        tsv=PRE_BOWTIE2 / "{genome}" / "{sample}.{library}.stats.tsv",
    log:
        PRE_BOWTIE2 / "{genome}" / "{sample}.{library}.stats.log",
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
