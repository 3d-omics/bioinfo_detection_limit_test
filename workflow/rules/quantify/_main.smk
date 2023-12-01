include: "_functions.smk"
include: "samtools.smk"
include: "coverm.smk"
include: "bowtie2.smk"


rule quantify:
    """Run the stats analyses: nonpareil and coverm"""
    input:
        rules.quantify__coverm.input,
