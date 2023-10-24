include: "functions.smk"
include: "fastp.smk"
include: "bowtie2_hosts.smk"
include: "bowtie2_mags.smk"


rule pre:
    input:
        rules.fastp.input,
        rules.bowtie2_hosts.input,
        rules.bowtie2_mags.input,
