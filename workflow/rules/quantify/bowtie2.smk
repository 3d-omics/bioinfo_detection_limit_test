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
        "_env.yml"
    params:
        extra=params["quantify"]["bowtie2"]["extra"],
    threads: 24
    resources:
        mem_mb=double_ram(32),
        runtime=24 * 60,
    retries: 5
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.reference} \
            {output.mock} \
        2> {log} 1>&2
        """


rule _quantify__bowtie2__map:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=get_host_clean_forward,
        reverse_=get_host_clean_reverse,
        mock=QUANT_BOWTIE2 / "{mag_catalogue}_index",
        reference=REFERENCE / "mags" / "{mag_catalogue}.fa.gz",
    output:
        cram=QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.cram",
        crai=QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.cram.crai",
    log:
        QUANT_BOWTIE2 / "{mag_catalogue}" / "{sample_id}.{library_id}.log",
    conda:
        "_env.yml"
    threads: 24
    resources:
        mem_mb=double_ram(32),
        runtime=24 * 60,
    retries: 5
    params:
        extra=params["quantify"]["bowtie2"]["extra"],
        samtools_mem=params["quantify"]["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    shell:
        """
        find \
            $(dirname {output.cram}) \
            -name "$(basename {output.cram}).tmp.*.bam" \
            -delete \
        2> {log} 1>&2

        ( bowtie2 \
            -x {input.mock} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
            --rg '{params.rg_extra}' \
            --rg-id '{params.rg_id}' \
            --threads {threads} \
            {params.extra} \
        | samtools sort \
            --reference {input.reference} \
            --threads {threads} \
            --write-index \
            -M \
            -T {output.cram} \
            -l 9 \
            -m {params.samtools_mem} \
            -o {output.cram} \
        ) 2>> {log} 1>&2
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
