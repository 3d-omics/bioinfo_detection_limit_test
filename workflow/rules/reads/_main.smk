include: "_functions.smk"


rule _reads__link:
    """Make a link to the original forward file, with a prettier name than default"""
    input:
        forward_=get_forward,
        reverse_=get_reverse,
    output:
        forward_=temp(READS / "{sample}.{library}_1.fq.gz"),
        reverse_=temp(READS / "{sample}.{library}_2.fq.gz"),
    log:
        READS / "{sample}.{library}.log",
    conda:
        "_env.yml"
    shell:
        """
        ln --symbolic $(readlink --canonicalize {input.forward_}) {output.forward_} 2> {log}
        ln --symbolic $(readlink --canonicalize {input.reverse_}) {output.reverse_} 2> {log}
        """


rule reads__link:
    """Link all reads in the samples.tsv"""
    input:
        [
            READS / f"{sample}.{library}_{end}.fq.gz"
            for sample, library in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],


rule reads__fastqc:
    """Run fastqc on all raw reads"""
    input:
        [
            READS / f"{sample}.{library}_{end}_fastqc.{extension}"
            for sample, library in SAMPLE_LIBRARY
            for end in ["1", "2"]
            for extension in ["html", "zip"]
        ],


rule reads:
    """Link all reads and run fastqc on them"""
    input:
        rules.reads__link.input,
        rules.reads__fastqc.input,


localrules:
    _reads__link,
