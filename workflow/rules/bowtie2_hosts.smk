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


rule bowtie2_hosts_map_pe_one:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=get_input_forward_for_host_mapping,
        reverse_=get_input_reverse_for_host_mapping,
        unpaired=get_input_single_for_host_mapping,
        mock=BOWTIE2_HOSTS / "{genome}_index",
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        cram=BOWTIE2_HOSTS / "{genome}/{sample}.{library}_pe.cram",
        crai=BOWTIE2_HOSTS / "{genome}/{sample}.{library}_pe.cram.crai",
    log:
        BOWTIE2_HOSTS / "{genome}/{sample}.{library}_pe.log",
    params:
        extra=params["bowtie2"]["extra"],
        samtools_mem=params["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    threads: 24
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
            -U {input.unpaired} \
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
            {params.extra} \
        | samtools sort \
            --threads {threads} \
            -m {params.samtools_mem} \
        | samtools rmdup - - \
        | samtools view \
            --reference {input.reference} \
            --output {output.cram} \
            --output-fmt cram,level=9,nthreads={threads} \
            --write-index \
        ) 2> {log} 1>&2
        """


rule bowtie2_hosts_map_se_one:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        single=get_input_single_for_host_mapping,
        mock=BOWTIE2_HOSTS / "{genome}_index",
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        cram=BOWTIE2_HOSTS / "{genome}/{sample}.{library}_se.cram",
        crai=BOWTIE2_HOSTS / "{genome}/{sample}.{library}_se.cram.crai",
    log:
        BOWTIE2_HOSTS / "{genome}/{sample}.{library}_se.log",
    params:
        extra=params["bowtie2"]["extra"],
        samtools_mem=params["bowtie2"]["samtools"]["mem_per_thread"],
        rg_id=compose_rg_id,
        rg_extra=compose_rg_extra,
    threads: 24
    conda:
        "../envs/bowtie2.yml"
    resources:
        mem_mb=32 * 1024,
        runtime=24 * 60,
    shell:
        """
        (bowtie2 \
            -x {input.mock} \
            -U {input.single} \
            --threads {threads} \
            --rg-id '{params.rg_id}' \
            --rg '{params.rg_extra}' \
            {params.extra} \
        | samtools sort \
            --threads {threads} \
            -m {params.samtools_mem} \
        | samtools rmdup - - \
        | samtools view \
            --reference {input.reference} \
            --output {output.cram} \
            --output-fmt cram,level=9,nthreads={threads} \
            --write-index \
        ) 2> {log} 1>&2
        """


rule bowtie2_hosts_extract_pe_one:
    """
    Keep only pairs unmapped to the human reference genome, sort by name rather
    than by coordinate, and convert to FASTQ.
    """
    input:
        cram=BOWTIE2_HOSTS / "{genome}/{sample}.{library}_pe.cram",
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        forward_=BOWTIE2_HOSTS / "non{genome}/{sample}.{library}_1.fq.gz",
        reverse_=BOWTIE2_HOSTS / "non{genome}/{sample}.{library}_2.fq.gz",
    log:
        BOWTIE2_HOSTS / "non{genome}/{sample}.{library}.log",
    conda:
        "../envs/bowtie2.yml"
    threads: 24
    resources:
        runtime=1 * 60,
        mem_mb=32 * 1024,
    shell:
        """
        (samtools view \
            --reference {input.reference} \
            --threads {threads} \
            -u \
            -o /dev/stdout \
            -f 12 \
            {input.cram} \
        | samtools fastq \
            -1 {output.forward_} \
            -2 {output.reverse_} \
            -0 /dev/null \
            -c 9 \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule bowtie2_hosts_extract_se_one:
    """
    Keep only pairs unmapped to the human reference genome, sort by name rather
    than by coordinate, and convert to FASTQ.
    """
    input:
        cram=BOWTIE2_HOSTS / "{genome}/{sample}.{library}_se.cram",
        reference=REFERENCE / "{genome}.fa.gz",
    output:
        forward_=BOWTIE2_HOSTS / "non{genome}/{sample}.{library}_se.fq.gz",
    log:
        BOWTIE2_HOSTS / "non{genome}/{sample}.{library}.log",
    conda:
        "../envs/bowtie2.yml"
    threads: 24
    resources:
        runtime=1 * 60,
        mem_mb=32 * 1024,
    shell:
        """
        (samtools view \
            --reference {input.reference} \
            --threads {threads} \
            -u \
            -o /dev/stdout \
            -f 12 \
            {input.cram} \
        | samtools fastq \
            -1 {output.forward_} \
            -2 /dev/null \
            -0 /dev/null \
            -c 9 \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule bowtie2_hosts_extract:
    """Run bowtie2_extract_nonchicken_one for all libraries"""
    input:
        [
            BOWTIE2_HOSTS / f"non{genome}/{sample}.{library}_{end}.fq.gz"
            for genome in [LAST_HOST]
            for sample, library in SAMPLE_LIB_PE
            for end in ["1", "2"]
        ],
        [
            BOWTIE2_HOSTS / f"non{genome}/{sample}.{library}_{end}.fq.gz"
            for genome in [LAST_HOST]
            for sample, library in SAMPLE_LIB_SE
            for end in ["se"]
        ],


rule bowtie2_hosts_report_all:
    """Generate bowtie2 reports for all libraries:
    - samtools stats
    - samtools flagstats
    - samtools idxstats
    """
    input:
        [
            BOWTIE2_HOSTS / f"{genome}/{sample}.{library}_{library_type}.{report}"
            for genome in HOST_NAMES
            for sample, library in SAMPLE_LIB_PE
            for library_type in ["pe", "se"]
            for report in BAM_REPORTS
        ],


rule bowtie2_host:
    """Run bowtie2 on all libraries and generate reports"""
    input:
        rules.bowtie2_hosts_report_all.input,
        rules.bowtie2_hosts_extract.input,
