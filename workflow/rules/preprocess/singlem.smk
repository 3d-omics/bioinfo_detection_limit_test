rule _preprocess__singlem__pipe:
    """Run singlem over one sample

    Note: SingleM asks in the documentation for the raw reads. Here we are
    passing it the non-host and trimmed ones.
    """
    input:
        forward_=get_host_clean_forward,
        reverse_=get_host_clean_reverse,
        data=features["databases"]["singlem"],
    output:
        archive_otu_table=SINGLEM / "pipe" / "{sample}.{library}.archive.json",
        otu_table=SINGLEM / "pipe" / "{sample}.{library}.otu_table.tsv",
        condense=SINGLEM / "pipe" / "{sample}.{library}.condense.tsv",
    log:
        SINGLEM / "pipe" / "{sample}.{library}.log",
    conda:
        "_env.yml"
    resources:
        runtime=4 * 60,
        mem_mb=8 * 1024,
    params:
        input_string=get_input_string_for_stats_singlem_pipe_one,
    shell:
        """
        singlem pipe \
            {params.input_string} \
            --otu-table {output.otu_table} \
            --archive-otu-table {output.archive_otu_table} \
            --taxonomic-profile {output.condense} \
            --metapackage {input.data} \
            --assignment-threads {threads} \
        2> {log} 1>&2
        """


rule _preprocess__singlem__condense:
    """Aggregate all the singlem results into a single table"""
    input:
        archive_otu_tables=[
            SINGLEM / "pipe" / f"{sample}.{library}.archive.json"
            for sample, library in SAMPLE_LIB
        ],
        data=features["databases"]["singlem"],
    output:
        condense=SINGLEM / "singlem.tsv",
    log:
        SINGLEM / "singlem.log",
    conda:
        "_env.yml"
    params:
        input_dir=SINGLEM,
    shell:
        """
        singlem condense \
            --input-archive-otu-tables {input.archive_otu_tables} \
            --taxonomic-profile {output.condense} \
            --metapackage {input.data} \
        2> {log} 1>&2
        """


rule preprocess__singlem:
    """Run all stats singlem rules"""
    input:
        SINGLEM / "singlem.tsv",
