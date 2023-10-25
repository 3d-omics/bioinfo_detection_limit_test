rule bowtie2_mags_build:
    """Build bowtie2 index for the mag reference

    Let the script decide to use a small or a large index based on the size of
    the reference genome.
    """
    input:
        reference=REFERENCE / "mags/{mag_catalogue}.fa.gz",
    output:
        mock=touch(BOWTIE2_MAGS / "{mag_catalogue}_index"),
    log:
        BOWTIE2_MAGS / "{mag_catalogue}_index.log",
    conda:
        "pre.yml"
    params:
        extra=params["pre"]["bowtie2"]["extra"],
    threads: 8
    resources:
        mem_mb=32 * 1024,
        runtime=24 * 60,
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {params.extra} \
            {input.reference} \
            {output.mock} \
        2> {log} 1>&2
        """


rule bowtie2_mags_map_one_library_to_one_catalogue:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=get_input_forward_for_mag_mapping,
        reverse_=get_input_reverse_for_mag_mapping,
        mock=BOWTIE2_MAGS / "{mag_catalogue}_index",
        reference=REFERENCE / "mags/{mag_catalogue}.fa.gz",
    output:
        cram=BOWTIE2_MAGS / "{mag_catalogue}/{sample}.{library}.cram",
        crai=BOWTIE2_MAGS / "{mag_catalogue}/{sample}.{library}.cram.crai",
    log:
        BOWTIE2_MAGS / "{mag_catalogue}/{sample}.{library}.log",
    conda:
        "pre.yml"
    threads: 4
    resources:
        mem_mb=32 * 1024,
        runtime=24 * 60,
    params:
        is_paired=is_paired,
        extra=params["pre"]["bowtie2"]["extra"],
        samtools_mem=params["pre"]["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
        input_string=compose_input_string_for_bowtie2_mags_map_one_library_to_one_catalogue,
    shell:
        """
        ( bowtie2 \
            -x {input.mock} \
            {params.input_string} \
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
            {params.extra} \
        | samtools sort \
            -l 9 \
            -M \
            -m {params.samtools_mem} \
            -o {output.cram} \
            --reference {input.reference} \
            --threads {threads} \
            --write-index \
        ) 2> {log} 1>&2
        """


rule bowtie2_mags_map_all_libraries_to_all_mags:
    """Run bowtie2_map_mags_one for all PE libraries"""
    input:
        [
            BOWTIE2_MAGS / f"{mag_catalogue}/{sample}.{library}.cram"
            for sample, library in SAMPLE_LIB
            for mag_catalogue in MAG_CATALOGUES
        ],


rule bowtie2_mags_report_all:
    """Generate bowtie2 reports for all PE libraries:
    - samtools stats
    - samtools flagstats
    - samtools idxstats
    """
    input:
        [
            BOWTIE2_MAGS / f"{mag_catalogue}/{sample}.{library}.{report}"
            for sample, library in SAMPLE_LIB
            for report in BAM_REPORTS
            for mag_catalogue in MAG_CATALOGUES
        ],


rule bowtie2_mags:
    """Run bowtie2 on all libraries and generate reports"""
    input:
        rules.bowtie2_mags_report_all.input,
        rules.bowtie2_mags_map_all_libraries_to_all_mags.input,
