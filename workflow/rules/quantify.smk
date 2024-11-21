include: "quantify/mags.smk"
include: "quantify/bowtie2.smk"
include: "quantify/coverm.smk"
include: "quantify/multiqc.smk"


rule quantify__all:
    input:
        rules.quantify__mags__all.input,
        rules.quantify__bowtie2__all.input,
        rules.quantify__coverm__all.input,
        rules.quantify__multiqc__all.input,
