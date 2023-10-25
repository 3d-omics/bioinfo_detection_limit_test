rule bowtie2_hosts_build:
    """Build bowtie2 index for the human reference

    Let the script decide to use a small or a large index based on the size of
    the reference genome.
    """
    input:
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        mock=touch(BOWTIE2_HOSTS / "{genome}_index"),
    log:
        BOWTIE2_HOSTS / "{genome}_index.log",
    conda:
        "pre.yml"
    params:
        extra=params["pre"]["bowtie2"]["extra"],
    threads: 8
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


rule bowtie2_hosts_map_one:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=get_input_forward_for_host_mapping,
        reverse_=get_input_reverse_for_host_mapping,
        mock=BOWTIE2_HOSTS / "{genome}_index",
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        cram=BOWTIE2_HOSTS / "{genome}" / "{sample}.{library}.cram",
        crai=BOWTIE2_HOSTS / "{genome}" / "{sample}.{library}.cram.crai",
    log:
        BOWTIE2_HOSTS / "{genome}" / "{sample}.{library}_pe.log",
    params:
        extra=params["pre"]["bowtie2"]["extra"],
        samtools_mem=params["pre"]["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
        input_string=compose_input_string_for_bowtie2_hosts_map_one,
        rmdup_string=compose_rmdup_string_for_bowtie2_hosts_map_one,
    threads: 24
    conda:
        "pre.yml"
    resources:
        mem_mb=double_ram(32),
        runtime=24 * 60,
    retries: 5
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
            --threads {threads} \
            -m {params.samtools_mem} \
        | samtools rmdup \
            {params.rmdup_string} \
            - - \
        | samtools view \
            --reference {input.reference} \
            --output {output.cram} \
            --output-fmt cram,level=9,nthreads={threads} \
            --write-index \
        ) 2> {log} 1>&2
        """


rule bowtie2_hosts_extract_one:
    """
    Keep only pairs unmapped to the human reference genome, sort by name rather
    than by coordinate, and convert to FASTQ.
    """
    input:
        cram=BOWTIE2_HOSTS / "{genome}" / "{sample}.{library}.cram",
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        forward_=BOWTIE2_HOSTS / "non{genome}" / "{sample}.{library}_1.fq.gz",
        reverse_=touch(BOWTIE2_HOSTS / "non{genome}" / "{sample}.{library}_2.fq.gz"),
    log:
        BOWTIE2_HOSTS / "non{genome}" / "{sample}.{library}.log",
    conda:
        "pre.yml"
    threads: 24
    resources:
        runtime=1 * 60,
        mem_mb=double_ram(32),
    params:
        output_string=compose_output_string_for_bowtie2_hosts_extract_one,
        filter_int=compose_filter_int_for_bowtie2_hosts_extract_one,
    retries: 5
    shell:
        """
        ( samtools view \
            --reference {input.reference} \
            --threads {threads} \
            -u \
            -o /dev/stdout \
            -f {params.filter_int} \
            {input.cram} \
        | samtools collate \
            -O \
            -u \
            -f \
            --reference {input.reference} \
            -@ {threads} \
            - \
        | samtools fastq \
            {params.output_string} \
            -c 9 \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule bowtie2_hosts_extract:
    """Run bowtie2_extract_nonchicken_one for all PE libraries"""
    input:
        [
            BOWTIE2_HOSTS / f"non{genome}" / f"{sample}.{library}_{end}.fq.gz"
            for genome in [LAST_HOST]
            for sample, library in SAMPLE_LIB_PE
            for end in ["1", "2"]
        ],
        [
            BOWTIE2_HOSTS / f"non{genome}" / f"{sample}.{library}_{end}.fq.gz"
            for genome in [LAST_HOST]
            for sample, library in SAMPLE_LIB_SE
            for end in ["1"]
        ],


rule bowtie2_hosts_report:
    """Generate bowtie2 reports for all libraries:
    - samtools stats
    - samtools flagstats
    - samtools idxstats
    """
    input:
        [
            BOWTIE2_HOSTS / genome / f"{sample}.{library}.{report}"
            for genome in HOST_NAMES
            for sample, library in SAMPLE_LIB
            for report in BAM_REPORTS
        ],


rule bowtie2_hosts:
    """Run bowtie2 on all libraries and generate reports"""
    input:
        rules.bowtie2_hosts_report.input,
        rules.bowtie2_hosts_extract.input,
