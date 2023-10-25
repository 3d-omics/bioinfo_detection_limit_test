

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
        "stats.yml"
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
        tsv=touch(
            STATS_COVERM / "{mag_catalogue}/genome/{method}/{sample}.{library}.tsv"
        ),
    conda:
        "stats.yml"
    log:
        STATS_COVERM / "{mag_catalogue}/genome/{method}/{sample}.{library}.log",
    params:
        method="{method}",
        min_covered_fraction=params["stats"]["coverm"]["genome"]["min_covered_fraction"],
        separator=params["stats"]["coverm"]["genome"]["separator"],
    shell:
        """
        coverm genome \
            --bam-files {input.bam} \
            --methods {params.method} \
            --separator {params.separator} \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output} 2> {log} || true \
        """


rule stats_coverm_genome_aggregate_one_mag_catalogue:
    """Aggregate all the nonpareil results into a single table"""
    input:
        get_coverm_genome_tsv_files_for_aggregation,
    output:
        STATS / "coverm_genome_{mag_catalogue}.{method}.tsv",
    log:
        STATS / "coverm_genome_{mag_catalogue}.{method}.log",
    conda:
        "stats.yml"
    params:
        input_dir=lambda wildcards: STATS_COVERM
        / wildcards.mag_catalogue
        / "genome"
        / wildcards.method,
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
            STATS / f"coverm_genome_{mag_catalogue}.{method}.tsv"
            for mag_catalogue in MAG_CATALOGUES
            for method in params["stats"]["coverm"]["genome"]["methods"]
        ],


rule stats_coverm_contig_one_library_one_mag_catalogue:
    """Run coverm contig for one library and one mag catalogue"""
    input:
        bam=STATS_COVERM / "{mag_catalogue}/bams/{sample}.{library}.bam",
    output:
        tsv=STATS_COVERM / "{mag_catalogue}/contig/{method}/{sample}.{library}.tsv",
    conda:
        "stats.yml"
    log:
        STATS_COVERM / "{mag_catalogue}/contig/{method}/{sample}.{library}.log",
    params:
        method=lambda wildcards: wildcards.method,
    shell:
        """
        coverm contig \
            --bam-files {input.bam} \
            --methods {params.method} \
            --proper-pairs-only \
        > {output} 2> {log} || true
        """


rule stats_coverm_contig_aggregate_mag_catalogue:
    """Aggregate all the nonpareil results into a single table"""
    input:
        get_coverm_contig_tsv_files_for_aggregation,
    output:
        STATS / "coverm_contig_{mag_catalogue}.{method}.tsv",
    log:
        STATS / "coverm_contig_{mag_catalogue}.{method}.log",
    conda:
        "stats.yml"
    params:
        input_dir=lambda wildcards: STATS_COVERM
        / wildcards.mag_catalogue
        / "contig"
        / wildcards.method,
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
            STATS / f"coverm_contig_{mag_catalogue}.{method}.tsv"
            for mag_catalogue in MAG_CATALOGUES
            for method in params["stats"]["coverm"]["contig"]["methods"]
        ],


rule stats_coverm:
    """Run both coverm overall and contig"""
    input:
        rules.stats_coverm_genome.input,
        rules.stats_coverm_contig.input,
