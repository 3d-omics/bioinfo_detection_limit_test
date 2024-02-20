rule helpers__samtools__crai:
    input:
        "{prefix}.cram",
    output:
        "{prefix}.cram.crai",
    log:
        "{prefix}.cram.crai.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools index {input} 2> {log} 1>&2"


rule helpers__samtools__flagstats_cram:
    """Compute flagstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        txt="{prefix}.flagstats.txt",
    log:
        "{prefix}.flagstats.log",
    conda:
        "__environment__.yml"
    resources:
        mem_mb=4 * 1024,
    shell:
        "samtools flagstats {input.cram} > {output.txt} 2> {log}"


rule helpers__samtools__idxstats_cram:
    """Compute idxstats for a cram"""
    input:
        cram="{prefix}.cram",
        crai="{prefix}.cram.crai",
    output:
        tsv="{prefix}.idxstats.tsv",
    log:
        "{prefix}.idxstats.log",
    conda:
        "__environment__.yml"
    resources:
        mem_mb=4 * 1024,
    shell:
        "samtools idxstats {input.cram} > {output.tsv} 2> {log}"


rule helpers__samtools__index_fa_gz:
    """Index a fa.gz file"""
    input:
        "{prefix}.fa.gz",
    output:
        "{prefix}.fa.gz.fai",
    log:
        "{prefix}.fa.gz.fai.log",
    conda:
        "__environment__.yml"
    shell:
        "samtools faidx {input} 2> {log} 1>&2"
