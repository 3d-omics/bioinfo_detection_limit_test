rule stats_nonpareil_one:
    """
    Note: Nonpareil only ask for one of the pair-end reads
    Note2: it has to be fastq. The process substitution trick does not work
    """
    input:
        forward_=BOWTIE2_NONCHICKEN / "{sample}.{library}_1.fq.gz",
    output:
        forward_fq=temp(STATS_NONPAREIL / "{sample}.{library}_1.fq"),
        npa=STATS_NONPAREIL / "{sample}.{library}.npa",
        npc=STATS_NONPAREIL / "{sample}.{library}.npc",
        npl=STATS_NONPAREIL / "{sample}.{library}.npl",
        npo=STATS_NONPAREIL / "{sample}.{library}.npo",
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
        2>> {log} 1>&2
        """


rule stats_nonpareil_all:
    input:
        [
            STATS_NONPAREIL / f"{sample}.{library}.{extension}"
            for extension in ["npa", "npc", "npl", "npo"]
            for sample, library in SAMPLE_LIB
        ],


rule stats_singlem_one:
    """
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
    input:
        [
            STATS_SINGLEM / f"{sample}.{library}.otu_table.tsv"
            for sample, library in SAMPLE_LIB
        ],


rule stats_cram_to_mapped_bam:
    input:
        cram=BOWTIE2_MAGS / "{sample}.{library}.cram",
        reference=REFERENCE / "mags.fa.gz",
    output:
        bam=temp(STATS_COVERM / "{sample}.{library}.bam"),
    log:
        STATS_COVERM / "{sample}.{library}.log",
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


rule stats_coverm_overall:
    input:
        bams=[
            STATS_COVERM / f"{sample}.{library}.bam" for sample, library in SAMPLE_LIB
        ],
    output:
        STATS / "coverm_overall.tsv",
    log:
        STATS / "coverm_overall.log",
    conda:
        "../envs/stats.yml"
    params:
        methods=params["coverm"]["genome"]["methods"],
        min_covered_fraction=params["coverm"]["genome"]["min_covered_fraction"],
        separator=params["coverm"]["genome"]["separator"],
    threads: 24
    resources:
        runtime=24 * 60,
        mem_mb=32 * 1024,
    shell:
        """
        coverm genome \
            --bam-files {input.bams} \
            --methods {params.methods} \
            --separator ^ \
            --threads {threads} \
            --min-covered-fraction {params.min_covered_fraction} \
        > {output} \
        2> {log}
        """


rule stats_coverm_contig:
    input:
        bams=[
            STATS_COVERM / f"{sample}.{library}.bam" for sample, library in SAMPLE_LIB
        ],
    output:
        STATS / "coverm_contig.tsv",
    log:
        STATS / "coverm_contig.log",
    conda:
        "../envs/stats.yml"
    params:
        methods=params["coverm"]["contig"]["methods"],
    threads: 24
    resources:
        runtime=24 * 60,
        mem_mb=32 * 1024,
    shell:
        """
        coverm contig \
            --threads {threads} \
            --bam-files {input.bams} \
            --methods {params.methods} \
            --proper-pairs-only \
        > {output} \
        2> {log}
        """


rule stats_coverm_all:
    input:
        rules.stats_coverm_overall.output,
        rules.stats_coverm_contig.output,


rule stats:
    input:
        rules.stats_nonpareil_all.input,
        rules.stats_singlem_all.input,
        rules.stats_coverm_all.input,
