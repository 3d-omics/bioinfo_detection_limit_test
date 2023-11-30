rule _stats__nonpareil__run:
    """Run nonpareil over one sample

    Note: Nonpareil only ask for one of the pair-end reads
    Note2: it has to be fastq. The process substitution trick does not work
    Note3: in case that nonpareil fails for low coverage samples, it creates
    empty files
    """
    input:
        forward_=get_input_forward_for_stats,
    output:
        forward_fq=temp(NONPAREIL / "{sample}.{library}_1.fq"),
        npa=touch(NONPAREIL / "{sample}.{library}.npa"),
        npc=touch(NONPAREIL / "{sample}.{library}.npc"),
        npl=touch(NONPAREIL / "{sample}.{library}.npl"),
        npo=touch(NONPAREIL / "{sample}.{library}.npo"),
    log:
        NONPAREIL / "{sample}.{library}.log",
    conda:
        "_env.yml"
    params:
        prefix=compose_prefix_for_nonpareil,
    resources:
        runtime=24 * 60,
    shell:
        """
        gzip -dc {input.forward_} > {output.forward_fq} 2> {log}

        nonpareil \
            -s {output.forward_fq} \
            -T kmer \
            -b {params.prefix} \
            -f fastq \
            -t {threads} \
        2>> {log} 1>&2 || true
        """


rule stats__nonpareil__run:
    """Run stats_nonpareil_one for all the samples"""
    input:
        [
            NONPAREIL / f"{sample}.{library}.{extension}"
            for extension in ["npa", "npc", "npl", "npo"]
            for sample, library in SAMPLE_LIB
        ],


rule _stats_nonpareil__aggregate:
    """Aggregate all the nonpareil results into a single table"""
    input:
        rules.stats__nonpareil__run.input,
    output:
        STATS / "nonpareil.tsv",
    log:
        STATS / "nonpareil.log",
    conda:
        "_env.yml"
    params:
        input_dir=NONPAREIL,
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_nonpareil.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule stats__nonpareil:
    input:
        rules._stats_nonpareil__aggregate.output,
