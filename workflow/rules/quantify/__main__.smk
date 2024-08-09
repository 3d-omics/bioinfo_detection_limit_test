include: "__functions__.smk"
include: "index.smk"
include: "bowtie2.smk"
include: "samtools.smk"
include: "coverm.smk"


rule quantify:
    """Quantify MAG abundances"""
    input:
        rules.quantify__coverm.input,
