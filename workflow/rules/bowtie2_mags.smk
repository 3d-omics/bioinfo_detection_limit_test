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
        "../envs/bowtie2.yml"
    params:
        extra=params["bowtie2"]["extra"],
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
        forward_=BOWTIE2_NONCHICKEN / "{sample}.{library}_1.fq.gz",
        reverse_=BOWTIE2_NONCHICKEN / "{sample}.{library}_2.fq.gz",
        mock=BOWTIE2_MAGS / "{mag_catalogue}_index",
        reference=REFERENCE / "mags/{mag_catalogue}.fa.gz",
    output:
        cram=BOWTIE2_MAGS / "{mag_catalogue}/{sample}.{library}.cram",
    log:
        BOWTIE2_MAGS / "{mag_catalogue}/{sample}.{library}.log",
    params:
        extra=params["bowtie2"]["extra"],
        samtools_mem=params["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    threads: 4
    conda:
        "../envs/bowtie2.yml"
    resources:
        mem_mb=32 * 1024,
        runtime=24 * 60,
    shell:
        """
        (bowtie2 \
            -x {input.mock} \
            -1 {input.forward_} \
            -2 {input.reverse_} \
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
        ) 2> {log} 1>&2
        """


rule bowtie2_mags_map_all_libraries_to_all_mags:
    """Run bowtie2_map_mags_one for all libraries"""
    input:
        [
            BOWTIE2_MAGS / f"{mag_catalogue}/{sample}.{library}.cram"
            for sample, library in SAMPLE_LIB
            for mag_catalogue in MAG_CATALOGUES
        ],


rule bowtie2_mags_report_all:
    """Generate bowtie2 reports for all libraries:
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
