rule _helpers__fastqc__run:
    """Run FastQC on a FASTQ file"""
    input:
        "{prefix}.fq.gz",
    output:
        html="{prefix}_fastqc.html",
        zip="{prefix}_fastqc.zip",
    conda:
        "__env__.yml"
    log:
        "{prefix}_fastqc.log",
    shell:
        "fastqc {input} 2> {log} 1>&2"
