rule _preprocess__nonpareil__run:
    """Run nonpareil over one sample

    Note: Nonpareil only ask for one of the pair-end reads
    Note2: it has to be fastq. The process substitution trick does not work
    Note3: in case that nonpareil fails for low coverage samples, it creates
    empty files
    """
    input:
        forward_=get_host_clean_forward,
    output:
        npa=touch(NONPAREIL / "run" / "{sample}.{library}.npa"),
        npc=touch(NONPAREIL / "run" / "{sample}.{library}.npc"),
        npl=touch(NONPAREIL / "run" / "{sample}.{library}.npl"),
        npo=touch(NONPAREIL / "run" / "{sample}.{library}.npo"),
    log:
        NONPAREIL / "run" / "{sample}.{library}.log",
    conda:
        "_env.yml"
    params:
        prefix=compose_prefix_for_nonpareil,
        forward_fq=lambda wildcards: NONPAREIL
        / "run"
        / "{wildcards.sample}.{wildcards.library}_1.fq",
    resources:
        runtime=24 * 60,
    shell:
        """
        gzip \
            --decompress \
            --stdout \
            {input.forward_} \
        > {params.forward_fq} 2> {log}

        nonpareil \
            -s {params.forward_fq} \
            -T kmer \
            -b {params.prefix} \
            -f fastq \
            -t {threads} \
        2>> {log} 1>&2 || true

        rm --force {params.forward_fq} 2>> {log} 1>&2
        """


rule preprocess__nonpareil__run:
    """Run stats_nonpareil_one for all the samples"""
    input:
        [
            NONPAREIL / "run" / f"{sample}.{library}.{extension}"
            for extension in ["npa", "npc", "npl", "npo"]
            for sample, library in SAMPLE_LIB
        ],


rule _preprocess__nonpareil__aggregate:
    """Aggregate all the nonpareil results into a single table"""
    input:
        rules.preprocess__nonpareil__run.input,
    output:
        NONPAREIL / "nonpareil.tsv",
    log:
        NONPAREIL / "nonpareil.log",
    conda:
        "_env.yml"
    params:
        input_dir=NONPAREIL / "run",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_nonpareil.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule preprocess__nonpareil:
    input:
        rules._preprocess__nonpareil__aggregate.output,
