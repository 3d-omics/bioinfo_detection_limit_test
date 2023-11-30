rule report_step_reads:
    """Collect all reports for the reads step"""
    input:
        rules.reads__fastqc.input,
    output:
        html=REPORT_STEP / "reads.html",
    log:
        REPORT_STEP / "reads.log",
    conda:
        "report.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --filename reads \
            --title reads \
            --force \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report_step_fastp:
    """Collect all reports for the fastp step"""
    input:
        rules.preprocess__fastp__report.input,
    output:
        html=REPORT_STEP / "fastp.html",
    log:
        REPORT_STEP / "fastp.log",
    conda:
        "report.yml"
    params:
        dir=REPORT_STEP,
    shell:
        """
        multiqc \
            --title fastp \
            --force \
            --filename fastp \
            --outdir {params.dir} \
            {input} \
        2> {log} 1>&2
        """


rule report_step_kraken2_one:
    """Collect all reports for the kraken2 step and ONE database"""
    input:
        get_report_step_kraken2_reports,
    output:
        html=REPORT_STEP / "kraken2_{kraken2_db}.html",
    log:
        REPORT_STEP / "kraken2_{kraken2_db}.log",
    conda:
        "report.yml"
    params:
        dir=REPORT_STEP,
        title="kraken2_{kraken2_db}",
    shell:
        """
        multiqc \
            --title {params.title} \
            --force \
            --filename {params.title} \
            --outdir {params.dir} \
            --module kraken \
            {input} \
        2> {log} 1>&2
        """


rule report_step_kraken2_all:
    """Collect all reports for the kraken2 step and ALL databases"""
    input:
        [REPORT_STEP / f"kraken2_{kraken2_db}.html" for kraken2_db in KRAKEN2_DBS],


rule report_step_bowtie2_host_one:
    """Collect all reports for the bowtie2 step and ONE host"""
    input:
        reports=[
            BOWTIE2_HOSTS / genome / f"{sample}.{library}.{report}"
            for genome in ["{genome}"]
            for sample, library in SAMPLE_LIB
            for report in BAM_REPORTS
        ],
    output:
        html=REPORT_STEP / "bowtie2_host_{genome}.html",
    log:
        REPORT_STEP / "bowtie2_host_{genome}.log",
    conda:
        "report.yml"
    params:
        dir=REPORT_STEP,
        title="bowtie2_host_{genome}",
    shell:
        """
        multiqc \
            --title {params.title} \
            --force \
            --filename {params.title} \
            --outdir {params.dir} \
            --dirs \
            --dirs-depth 1 \
            {input.reports} \
        2> {log} 1>&2
        """


rule report_step_bowtie2_host_all:
    """Collect all reports for the bowtie2 step and ALL hosts"""
    input:
        [REPORT_STEP / f"bowtie2_host_{genome}.html" for genome in HOST_NAMES],


rule report_step_bowtie2_mags_one:
    """Collect all reports for the bowtie2 step when mapping to a mag catalogue"""
    input:
        reports=[
            BOWTIE2_MAGS / mag_catalogue / f"{sample}.{library}.{report}"
            for mag_catalogue in ["{mag_catalogue}"]
            for sample, library in SAMPLE_LIB_PE
            for report in BAM_REPORTS
        ],
    output:
        html=REPORT_STEP / "bowtie2_mags_{mag_catalogue}.html",
    log:
        REPORT_STEP / "bowtie2_mags_{mag_catalogue}.log",
    conda:
        "report.yml"
    params:
        dir=REPORT_STEP,
        title="bowtie2_mags_{mag_catalogue}",
    shell:
        """
        multiqc \
            --title {params.title} \
            --force \
            --filename {params.title} \
            --outdir {params.dir} \
            --dirs \
            --dirs-depth 1 \
            {input.reports} \
        2> {log} 1>&2
        """


rule report_step_bowtie2_mags_all:
    """Collect all reports for the bowtie2 step and ALL mag catalogues"""
    input:
        [
            REPORT_STEP / f"bowtie2_mags_{mag_catalogue}.html"
            for mag_catalogue in MAG_CATALOGUES
        ],


rule report_step:
    """Collect all per step reports for the pipeline"""
    input:
        rules.report_step_reads.output,
        rules.report_step_fastp.output,
        rules.report_step_kraken2_all.input,  # input!
        rules.report_step_bowtie2_host_all.input,
        rules.report_step_bowtie2_mags_all.output,


localrules:
    report_step_reads,
    report_step_fastp,
    report_step_kraken2_one,
    report_step_bowtie2_host_one,
    report_step_bowtie2_mags_one,
