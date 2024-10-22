include: "mags.smk"
include: "bowtie2.smk"
include: "coverm.smk"
include: "multiqc.smk"


rule quantify__all:
    input:
        rules.quantify__mags__all.input,
        rules.quantify__bowtie2__all.input,
        rules.quantify__coverm__all.input,
        rules.quantify__multiqc__all.input,
