rule stats_nonpareil_one:
    """Run nonpareil over one sample

    Note: Nonpareil only ask for one of the pair-end reads
    Note2: it has to be fastq. The process substitution trick does not work
    Note3: in case that nonpareil fails for low coverage samples, it creates
    empty files
    """
    input:
        forward_=BOWTIE2_NONCHICKEN / "{sample}.{library}_1.fq.gz",
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


rule stats_singlem_one:
    """Run singlem over one sample

    Note: SingleM asks in the documentation for the raw reads. Here we are
    passing it the non-host and trimmed ones.
    """
    input:
        forward_=BOWTIE2_NONCHICKEN / "{sample}.{library}_1.fq.gz",
        reverse_=BOWTIE2_NONCHICKEN / "{sample}.{library}_2.fq.gz",
    output:
        otu_table=STATS_SINGLEM / "{sample}.{library}.otu_table.tsv",
    log:
        STATS_SINGLEM / "{sample}.{library}.log",
    conda:
        "../envs/stats.yml"
    threads: params["singlem"]["threads"]
    resources:
        runtime=24 * 60,
        mem_mb=32 * 1024,
    shell:
        """
        singlem pipe \
            --forward {input.forward_} \
            --reverse {input.reverse_} \
            --otu_table {output.otu_table} \
            --threads {threads} \
        2> {log} 1>&2
        """


rule stats_singlem_all:
    """Run stats_singlem_one for all the samples"""
    input:
        [
            STATS_SINGLEM / f"{sample}.{library}.otu_table.tsv"
            for sample, library in SAMPLE_LIB
        ],


rule stats_singlem:
    """Aggregate all the singlem results into a single table"""
    input:
        rules.stats_singlem_all.input,
    output:
        STATS / "singlem.tsv",
    log:
        STATS / "singlem.log",
    conda:
        "../envs/stats_r.yml"
    params:
        input_dir=STATS_SINGLEM,
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_singlem.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log}
        """


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
    threads: 8
    resources:
        runtime=24 * 60,
        mem_mb=8 * 1024,
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
    input:
        bam=STATS_COVERM / "{mag_catalogue}/bams/{sample}.{library}.bam",
    output:
        tsv=STATS_COVERM / "{mag_catalogue}/genome/{sample}.{library}.tsv",
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
        > {output} \
        2> {log}
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
    input:
        [
            STATS / f"coverm_genome_{mag_catalogue}.tsv"
            for mag_catalogue in MAG_CATALOGUES
        ],


rule stats_coverm_contig_one_library_one_mag_catalogue:
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
        > {output} \
        2> {log}
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
    """Run all the stats rules: nonpareil, singlem, and coverm"""
    input:
        rules.stats_nonpareil.output,
        rules.stats_singlem.output,
        rules.stats_coverm.input,


localrules:
    stats_nonpareil,
    stats_singlem,
