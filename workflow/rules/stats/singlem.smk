rule stats_singlem_pipe_one:
    """Run singlem over one sample

    Note: SingleM asks in the documentation for the raw reads. Here we are
    passing it the non-host and trimmed ones.
    """
    input:
        forward_=get_input_forward_for_stats,
        reverse_=get_input_reverse_for_stats,
        data=features["singlem_database"],
    output:
        archive_otu_table=STATS_SINGLEM / "{sample}.{library}.archive.json",
        otu_table=STATS_SINGLEM / "{sample}.{library}.otu_table.tsv",
        condense=STATS_SINGLEM / "{sample}.{library}.condense.tsv",
    log:
        STATS_SINGLEM / "{sample}.{library}.log",
    conda:
        "stats.yml"
    threads: 4
    resources:
        runtime=4 * 60,
    params:
        is_paired=is_paired,
    shell:
        """
        if [[ {params.is_paired} = "True" ]] ; then
            singlem pipe \
                --forward {input.forward_} \
                --reverse {input.reverse_} \
                --otu-table {output.otu_table} \
                --archive-otu-table {output.archive_otu_table} \
                --taxonomic-profile {output.condense} \
                --metapackage {input.data} \
                --threads {threads} \
                --assignment-threads {threads} \
            2> {log} 1>&2 || true
        else
            singlem pipe \
                --forward {input.forward_} \
                --otu-table {output.otu_table} \
                --archive-otu-table {output.archive_otu_table} \
                --taxonomic-profile {output.condense} \
                --metapackage {input.data} \
                --threads {threads} \
                --assignment-threads {threads} \
            2> {log} 1>&2 || true
        fi
        """


rule stats_singlem_condense:
    """Aggregate all the singlem results into a single table"""
    input:
        archive_otu_tables=[
            STATS_SINGLEM / f"{sample}.{library}.archive.json"
            for sample, library in SAMPLE_LIB
        ],
        data=features["singlem_database"],
    output:
        condense=STATS / "singlem.tsv",
    log:
        STATS / "singlem.log",
    conda:
        "stats.yml"
    params:
        input_dir=STATS_SINGLEM,
    shell:
        """
        singlem condense \
            --input-archive-otu-tables {input.archive_otu_tables} \
            --taxonomic-profile {output.condense} \
            --metapackage {input.data} \
        2> {log} 1>&2
        """


rule stats_singlem:
    """Run all stats singlem rules"""
    input:
        STATS / "singlem.tsv",
