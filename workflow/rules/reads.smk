rule reads_link_forward:
    """Make a link to the original forward file, with a prettier name than default"""
    input:
        forward_=get_forward,
    output:
        forward_=temp(READS / "{sample}.{library}_1.fq.gz"),
    log:
        READS / "{sample}.{library}_1.log",
    conda:
        "../envs/reads.yml"
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.forward_}) {output.forward_}
        """


rule reads_link_reverse:
    """Make a link to the original reverse file, with a prettier name than default"""
    input:
        reverse_=get_reverse,
    output:
        reverse_=temp(READS / "{sample}.{library}_2.fq.gz"),
    log:
        READS / "{sample}.{library}_2.log",
    conda:
        "../envs/reads.yml"
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.reverse_}) {output.reverse_}
        """


rule reads_link_se:
    """Make a link to the original single end file, with a prettier name than default"""
    input:
        reverse_=get_forward,
    output:
        reverse_=temp(READS / "{sample}.{library}_se.fq.gz"),
    log:
        READS / "{sample}.{library}_se.log",
    conda:
        "../envs/reads.yml"
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.reverse_}) {output.reverse_}
        """


rule reads_link_all:
    """Link all reads in the samples.tsv"""
    input:
        [
            READS / f"{sample}.{library}_{end}.fq.gz"
            for sample, library in SAMPLE_LIB_PE
            for end in ["1", "2"]
        ],
        [READS / f"{sample}.{library}_se.fq.gz" for sample, library in SAMPLE_LIB_SE],


rule reads_fastqc_all:
    """Run fastqc on all raw reads"""
    input:
        [
            READS / f"{sample}.{library}_{end}_fastqc.{extension}"
            for sample, library in SAMPLE_LIB_PE
            for end in ["1", "2"]
            for extension in ["html", "zip"]
        ],
        [
            READS / f"{sample}.{library}_se_fastqc.{extension}"
            for sample, library in SAMPLE_LIB_SE
            for extension in ["html", "zip"]
        ],


rule reads:
    """Link all reads and run fastqc on them"""
    input:
        rules.reads_link_all.input,
        rules.reads_fastqc_all.input,


localrules:
    reads_link_forward,
    reads_link_reverse,
