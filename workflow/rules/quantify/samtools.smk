rule _quantify__samtools__stats__cram_mag:
    """Compute stats for a MAG cram"""
    input:
        cram=QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample}.{library}.cram",
        crai=QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample}.{library}.cram.crai",
        reference=REFERENCE / "mags" / "{mag_catalogue}.fa.gz",
    output:
        tsv=QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample}.{library}.stats.tsv",
    log:
        QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample}.{library}.stats.log",
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
