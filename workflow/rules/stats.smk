rule stats_nonpareil_one:
    """Run nonpareil over one sample

    Note: Nonpareil only ask for one of the pair-end reads
    Note2: it has to be fastq. The process substitution trick does not work
    Note3: in case that nonpareil fails for low coverage samples, it creates
    empty files
    """
    input:
        forward_=get_input_forward_for_stats,
    output:
        forward_fq=temp(STATS_NONPAREIL / "{sample}.{library}_1.fq"),
        npa=touch(STATS_NONPAREIL / "{sample}.{library}.npa"),
        npc=touch(STATS_NONPAREIL / "{sample}.{library}.npc"),
        npl=touch(STATS_NONPAREIL / "{sample}.{library}.npl"),
        npo=touch(STATS_NONPAREIL / "{sample}.{library}.npo"),
    log:
        STATS_NONPAREIL / "{sample}.{library}.log",
    conda:
        "../envs/stats.yml"
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


rule stats_nonpareil_all:
    """Run stats_nonpareil_one for all the samples"""
    input:
        [
            STATS_NONPAREIL / f"{sample}.{library}.{extension}"
            for extension in ["npa", "npc", "npl", "npo"]
            for sample, library in SAMPLE_LIB
        ],


rule stats_nonpareil:
    """Aggregate all the nonpareil results into a single table"""
    input:
        rules.stats_nonpareil_all.input,
    output:
        STATS / "nonpareil.tsv",
    log:
        STATS / "nonpareil.log",
    conda:
        "../envs/stats_r.yml"
    params:
        input_dir=STATS_NONPAREIL,
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_nonpareil.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


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
        "../envs/stats.yml"
    threads: 4
    resources:
        runtime=4 * 60,
    shell:
        """
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
        """


rule stats_singlem_pipe_all:
    """Run stats_singlem_one for all the samples"""
    input:
        [
            STATS_SINGLEM / f"{sample}.{library}.otu_table.tsv"
            for sample, library in SAMPLE_LIB
        ],


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
        "../envs/stats.yml"
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


rule stats_cram_to_mapped_bam:
    """Convert cram to bam

    Note: this step is needed because coverm probably does not support cram. The
    log from coverm shows failures to get the reference online, but nonetheless
    it works.
    """
    input:
        cram=BOWTIE2_MAGS / "{mag_catalogue}/{sample}.{library}.cram",
        reference=REFERENCE / "mags/{mag_catalogue}.fa.gz",
    output:
        bam=temp(STATS_COVERM / "{mag_catalogue}/bams/{sample}.{library}.bam"),
    log:
        STATS_COVERM / "{mag_catalogue}/bams/{sample}.{library}.log",
    conda:
        "../envs/samtools.yml"
    threads: 24
    resources:
        runtime=1 * 60,
        mem_mb=4 * 1024,
    shell:
        """
        samtools view \
            -F 4 \
            --threads {threads} \
            --reference {input.reference} \
            --output {output.bam} \
            --fast \
            {input.cram} \
        2> {log}
        """


rule stats_coverm_genome_one_library_one_mag_catalogue:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=STATS_COVERM / "{mag_catalogue}/bams/{sample}.{library}.bam",
    output:
        tsv=touch(STATS_COVERM / "{mag_catalogue}/genome/{sample}.{library}.tsv"),
    conda:
        "../envs/stats.yml"
    log:
        STATS_COVERM / "{mag_catalogue}/genome/{sample}.{library}.log",
    params:
        methods=params["coverm"]["genome"]["methods"],
        min_covered_fraction=params["coverm"]["genome"]["min_covered_fraction"],
        separator=params["coverm"]["genome"]["separator"],
    shell:
        """
        coverm genome \
            --bam-files {input.bam} \
            --methods {params.methods} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output} 2> {log} || true \
        """


rule stats_coverm_genome_aggregate_one_mag_catalogue:
    """Aggregate all the nonpareil results into a single table"""
    input:
        get_coverm_genome_tsv_files_for_aggregation,
    output:
        STATS / "coverm_genome_{mag_catalogue}.tsv",
    log:
        STATS / "coverm_genome_{mag_catalogue}.log",
    conda:
        "../envs/stats_r.yml"
    params:
        input_dir=lambda wildcards: STATS_COVERM / wildcards.mag_catalogue / "genome",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule stats_coverm_genome:
    """Run all rules to run coverm genome over all MAG catalogues"""
    input:
        [
            STATS / f"coverm_genome_{mag_catalogue}.tsv"
            for mag_catalogue in MAG_CATALOGUES
        ],


rule stats_coverm_contig_one_library_one_mag_catalogue:
    """Run coverm contig for one library and one mag catalogue"""
    input:
        bam=STATS_COVERM / "{mag_catalogue}/bams/{sample}.{library}.bam",
    output:
        tsv=STATS_COVERM / "{mag_catalogue}/contig/{sample}.{library}.tsv",
    conda:
        "../envs/stats.yml"
    log:
        STATS_COVERM / "{mag_catalogue}/contig/{sample}.{library}.log",
    params:
        methods=params["coverm"]["contig"]["methods"],
    shell:
        """
        coverm contig \
            --bam-files {input.bam} \
            --methods {params.methods} \
            --proper-pairs-only \
        > {output} 2> {log} || true
        """


rule stats_coverm_contig_aggregate_mag_catalogue:
    """Aggregate all the nonpareil results into a single table"""
    input:
        get_coverm_contig_tsv_files_for_aggregation,
    output:
        STATS / "coverm_contig_{mag_catalogue}.tsv",
    log:
        STATS / "coverm_contig_{mag_catalogue}.log",
    conda:
        "../envs/stats_r.yml"
    params:
        input_dir=lambda wildcards: STATS_COVERM / wildcards.mag_catalogue / "contig",
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule stats_coverm_contig:
    """Run all rules to run coverm contig over all MAG catalogues"""
    input:
        [
            STATS / f"coverm_contig_{mag_catalogue}.tsv"
            for mag_catalogue in MAG_CATALOGUES
        ],


rule stats_coverm:
    """Run both coverm overall and contig"""
    input:
        rules.stats_coverm_genome.input,
        rules.stats_coverm_contig.input,


rule stats:
    """Run the stats analyses: nonpareil and coverm"""
    input:
        rules.stats_nonpareil.output,
        # rules.stats_singlem.input,
        rules.stats_coverm.input,


rule stats_with_singlem:
    """Run the nonpareil, coverm and singlem"""
    input:
        rules.stats.input,
        rules.stats_singlem.input,


localrules:
    stats_nonpareil,
    stats_singlem,
