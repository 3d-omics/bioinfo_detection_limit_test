include: "_functions.smk"
include: "fastp.smk"
include: "bowtie2_hosts.smk"
include: "bowtie2_mags.smk"


rule preprocess:
    input:
        rules.preprocess__fastp.input,
        rules.preprocess__bowtie2__hosts.input,
        rules.preprocess__bowtie2__mags.input,
