rule _preprocess__bowtie2__build:
    """Build bowtie2 index for the human reference

    Let the script decide to use a small or a large index based on the size of
    the reference genome.
    """
    input:
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        mock=touch(PRE_BOWTIE2 / "{genome}_index"),
    log:
        PRE_BOWTIE2 / "{genome}_index.log",
    conda:
        "__environment__.yml"
    threads: 24
    resources:
        mem_mb=double_ram(params["preprocess"]["bowtie2"]["mem_gb"]),
        runtime=24 * 60,
    retries: 5
    cache: True
    shell:
        """
        bowtie2-build \
            --threads {threads} \
            {input.reference} \
            {output.mock} \
        2> {log} 1>&2
        """


rule _preprocess__bowtie2__map:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=get_input_forward_for_host_mapping,
        reverse_=get_input_reverse_for_host_mapping,
        mock=PRE_BOWTIE2 / "{genome}_index",
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        cram=PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.cram",
    log:
        PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}_pe.log",
    params:
        samtools_mem=params["preprocess"]["bowtie2"]["samtools_mem"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    threads: 24
    conda:
        "__environment__.yml"
    resources:
        mem_mb=double_ram(params["preprocess"]["bowtie2"]["mem_gb"]),
        runtime=24 * 60,
    retries: 5
    group:
        "preprocess"
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
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
        | samtools sort \
            --output-fmt CRAM \
            --reference {input.reference} \
            --threads {threads} \
            -T {output.cram} \
            -l 9 \
            -m {params.samtools_mem} \
            -o {output.cram}
        ) 2>> {log} 1>&2
        """


rule _preprocess__bowtie2__extract:
    """
    Keep only pairs unmapped to the human reference genome, sort by name rather
    than by coordinate, and convert to FASTQ.
    """
    input:
        cram=PRE_BOWTIE2 / "{genome}" / "{sample_id}.{library_id}.cram",
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        forward_=PRE_BOWTIE2 / "non{genome}" / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=touch(PRE_BOWTIE2 / "non{genome}" / "{sample_id}.{library_id}_2.fq.gz"),
    log:
        PRE_BOWTIE2 / "non{genome}" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    threads: 24
    resources:
        runtime=1 * 60,
        mem_mb=double_ram(params["preprocess"]["bowtie2"]["mem_gb"]),
    retries: 5
    group:
        "preprocess"
    shell:
        """
        ( samtools view \
            --reference {input.reference} \
            --threads {threads} \
            -u \
            -o /dev/stdout \
            -f 12 \
            {input.cram} \
        | samtools collate \
            -O \
            -u \
            -f \
            --reference {input.reference} \
            --threads {threads} \
            - \
        | samtools fastq \
            -1 {output.forward_} \
            -2 {output.reverse_} \
            -0 /dev/null \
            -c 9 \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule preprocess__bowtie2__extract:
    """Run bowtie2_extract_nonchicken_one for all PE libraries"""
    input:
        [
            PRE_BOWTIE2 / f"non{genome}" / f"{sample_id}.{library_id}_{end}.fq.gz"
            for genome in [LAST_HOST]
            for sample_id, library_id in SAMPLE_LIBRARY
            for end in ["1", "2"]
        ],


rule preprocess__bowtie2__report:
    """Generate bowtie2 reports for all libraries:
    - samtools stats
    - samtools flagstats
    - samtools idxstats
    """
    input:
        [
            PRE_BOWTIE2 / genome / f"{sample_id}.{library_id}.{report}"
            for genome in HOST_NAMES
            for sample_id, library_id in SAMPLE_LIBRARY
            for report in BAM_REPORTS
        ],


rule preprocess__bowtie2:
    """Run bowtie2 on all libraries and generate reports"""
    input:
        rules.preprocess__bowtie2__report.input,
        rules.preprocess__bowtie2__extract.input,
