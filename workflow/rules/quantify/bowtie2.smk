rule _quantify__bowtie2__build:
    """Build bowtie2 index for the mag reference

    Let the script decide to use a small or a large index based on the size of
    the reference genome.
    """
    input:
        reference=REFERENCE / "mags" / "{mag_catalogue}.fa.gz",
    output:
        mock=touch(QUANT_BOWTIE2 / "{mag_catalogue}_index"),
    log:
        QUANT_BOWTIE2 / "{mag_catalogue}_index.log",
    conda:
        "__environment__.yml"
    threads: 24
    resources:
        mem_mb=double_ram(params["preprocess"]["bowtie2"]["mem_gb"]),
        runtime=24 * 60,
    retries: 5
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.reference} \
            {output.mock} \
        2> {log} 1>&2
        """


rule _quantify__bowtie2__map:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        cram=get_host_clean_cram,
        mock=QUANT_BOWTIE2 / "{mag_catalogue}_index",
        reference=REFERENCE / "mags" / "{mag_catalogue}.fa.gz",
    output:
        cram=QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.cram",
    log:
        QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    threads: 24
    resources:
        mem_mb=double_ram(params["quantify"]["bowtie2"]["mem_gb"]),
        runtime=24 * 60,
    retries: 5
    params:
        samtools_mem=params["quantify"]["bowtie2"]["samtools_mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    group:
        "sample"
    cache: True
    shell:
        """
        ( samtools view \
            -u \
            -f 12 \
            {input.cram} \
        | samtools sort \
            -u \
            -n \
            --threads {threads} \
        | bowtie2 \
            -x {input.mock} \
            -b /dev/stdin \
            --align-paired-reads \
            --rg '{params.rg_extra}' \
            --rg-id '{params.rg_id}' \
            --threads {threads} \
        | samtools sort \
            --output-fmt CRAM \
            --reference {input.reference} \
            --threads {threads} \
            -T {output.cram} \
            -l 9 \
            -m {params.samtools_mem} \
            -o {output.cram} \
        ) 2> {log} 1>&2
        """


rule quantify__bowtie2__map:
    """Run bowtie2_map_mags_one for all PE libraries"""
    input:
        [
            QUANT_BOWTIE2 / mag_catalogue / f"{sample_id}.{library_id}.cram"
            for sample_id, library_id in SAMPLE_LIBRARY
            for mag_catalogue in MAG_CATALOGUES
        ],


rule quantify__bowtie2__report:
    """Generate bowtie2 reports for all PE libraries:
    - samtools stats
    - samtools flagstats
    - samtools idxstats
    """
    input:
        [
            QUANT_BOWTIE2 / mag_catalogue / f"{sample_id}.{library_id}.{report}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for report in BAM_REPORTS
            for mag_catalogue in MAG_CATALOGUES
        ],


rule quantify__bowtie2:
    """Run bowtie2 on all libraries and generate reports"""
    input:
        rules.quantify__bowtie2__report.input,
        rules.quantify__bowtie2__map.input,
