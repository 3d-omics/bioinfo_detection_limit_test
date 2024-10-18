include: "coverm_functions.smk"


rule quantify__coverm__genome:
    """Run coverm genome for one library and one mag catalogue"""
    input:
        bam=QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.bam",
    output:
        tsv=touch(
            COVERM / "{mag_catalogue}/genome/{method}/{sample_id}.{library_id}.tsv.gz"
        ),
    conda:
        "../environments/coverm.yml"
    log:
        COVERM / "{mag_catalogue}/genome/{method}/{sample_id}.{library_id}.log",
    params:
        method="{method}",
        min_covered_fraction=params["quantify"]["coverm"]["genome"][
            "min_covered_fraction"
        ],
        separator=params["quantify"]["coverm"]["genome"]["separator"],
    shell:
        """
        coverm genome \
            --bam-files {input.bam} \
            --methods {params.method} \
            --separator "{params.separator}" \
            --min-covered-fraction {params.min_covered_fraction} \
            --output-file {output.tsv} \
        2> {log} 1>&2
        """


rule quantify__coverm__genome__aggregate:
    """Aggregate all the nonpareil results into a single table"""
    input:
        get_coverm_genome_tsv_files_for_aggregation,
    output:
        COVERM / "coverm_genome_{mag_catalogue}.{method}.tsv.gz",
    log:
        COVERM / "coverm_genome_{mag_catalogue}.{method}.log",
    conda:
        "../environments/coverm.yml"
    params:
        input_dir=lambda w: COVERM / w.mag_catalogue / "genome" / w.method,
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule quantify__coverm__genome__all:
    """Run all rules to run coverm genome over all MAG catalogues"""
    input:
        [
            COVERM / f"coverm_genome_{mag_catalogue}.{method}.tsv.gz"
            for mag_catalogue in MAG_CATALOGUES
            for method in params["quantify"]["coverm"]["genome"]["methods"]
        ],


rule quantify__coverm__contig:
    """Run coverm contig for one library and one mag catalogue"""
    input:
        bam=QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.bam",
    output:
        tsv=COVERM / "{mag_catalogue}/contig/{method}/{sample_id}.{library_id}.tsv.gz",
    conda:
        "../environments/coverm.yml"
    log:
        COVERM / "{mag_catalogue}/contig/{method}/{sample_id}.{library_id}.log",
    params:
        method="{method}",
    shell:
        """
        coverm contig \
            --bam-files {input.bam} \
            --methods {params.method} \
            --proper-pairs-only \
            --output-file {output.tsv} \
        2> {log} 1>&2
        """


rule quantify__coverm__contig__aggregate:
    """Aggregate all the nonpareil results into a single table"""
    input:
        get_coverm_contig_tsv_files_for_aggregation,
    output:
        COVERM / "coverm_contig_{mag_catalogue}.{method}.tsv.gz",
    log:
        COVERM / "coverm_contig_{mag_catalogue}.{method}.log",
    conda:
        "../environments/coverm.yml"
    params:
        input_dir=lambda w: COVERM / w.mag_catalogue / "contig" / w.method,
    shell:
        """
        Rscript --no-init-file workflow/scripts/aggregate_coverm.R \
            --input-folder {params.input_dir} \
            --output-file {output} \
        2> {log} 1>&2
        """


rule quantify__coverm__contig__all:
    """Run all rules to run coverm contig over all MAG catalogues"""
    input:
        [
            COVERM / f"coverm_contig_{mag_catalogue}.{method}.tsv.gz"
            for mag_catalogue in MAG_CATALOGUES
            for method in params["quantify"]["coverm"]["contig"]["methods"]
        ],


rule quantify__coverm__all:
    """Run both coverm overall and contig"""
    input:
        rules.quantify__coverm__genome__all.input,
        rules.quantify__coverm__contig__all.input,
