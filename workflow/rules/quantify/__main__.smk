include: "__functions__.smk"
include: "index.smk"
include: "bowtie2.smk"
include: "samtools.smk"
include: "coverm.smk"


rule quantify:
    """Run the stats analyses: nonpareil and coverm"""
    input:
        rules.quantify__coverm.input,
