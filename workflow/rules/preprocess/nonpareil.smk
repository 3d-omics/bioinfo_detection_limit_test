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
        npa=touch(NONPAREIL / "run" / "{sample_id}.{library_id}.npa"),
        npc=touch(NONPAREIL / "run" / "{sample_id}.{library_id}.npc"),
        npl=touch(NONPAREIL / "run" / "{sample_id}.{library_id}.npl"),
        npo=touch(NONPAREIL / "run" / "{sample_id}.{library_id}.npo"),
    log:
        NONPAREIL / "run" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    params:
        prefix=compose_prefix_for_nonpareil,
        forward_fq=lambda w: NONPAREIL / "run" / f"{w.sample_id}.{w.library_id}_1.fq",
    resources:
        runtime=24 * 60,
    shell:
        """
        gzip \
            --decompress \
            --stdout \
            {input.forward_} \
        > {params.forward_fq} \
        2> {log}

        nonpareil \
            -s {params.forward_fq} \
            -T kmer \
            -b {params.prefix} \
            -f fastq \
            -t {threads} \
        2>> {log} \
        1>&2 || true

        rm --force --verbose {params.forward_fq} 2>> {log} 1>&2
        """


rule preprocess__nonpareil__run:
    """Run stats_nonpareil_one for all the samples"""
    input:
        [
            NONPAREIL / "run" / f"{sample_id}.{library_id}.{extension}"
            for extension in ["npa", "npc", "npl", "npo"]
            for sample_id, library_id in SAMPLE_LIBRARY
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
        "__environment__.yml"
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
