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


rule bowtie2_hosts_map_human_one:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=FASTP / "{sample}.{library}_1.fq.gz",
        reverse_=FASTP / "{sample}.{library}_2.fq.gz",
        mock=BOWTIE2_HOSTS / "human_index",
        reference=REFERENCE / "human.fa.gz",
    output:
        cram=BOWTIE2_HUMAN / "{sample}.{library}.cram",
    log:
        BOWTIE2_HUMAN / "{sample}.{library}.log",
    benchmark:
        BOWTIE2_HUMAN / "{sample}.{library}.bmk"
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


rule bowtie2_hosts_extract_nonhuman_one:
    """
    Keep only pairs unmapped to the human reference genome, sort by name rather
    than by coordinate, and convert to FASTQ.
    """
    input:
        cram=BOWTIE2_HUMAN / "{sample}.{library}.cram",
        reference=REFERENCE / "human.fa.gz",
    output:
        forward_=BOWTIE2_NONHUMAN / "{sample}.{library}_1.fq.gz",
        reverse_=BOWTIE2_NONHUMAN / "{sample}.{library}_2.fq.gz",
    log:
        BOWTIE2_NONHUMAN / "{sample}.{library}.log",
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
        | samtools sort \
            -n \
            -u \
            --threads {threads} \
        | samtools fastq \
            -1 {output.forward_} \
            -2 {output.reverse_} \
            -0 /dev/null \
            -c 9 \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule bowtie2_hosts_extract_nonhuman_all:
    """Run bowtie2_extract_nonhuman_one for all libraries"""
    input:
        [
            BOWTIE2_NONHUMAN / f"{sample}.{library}_{end}.fq.gz"
            for sample, library in SAMPLE_LIB
            for end in ["1", "2"]
        ],


rule bowtie2_hosts_map_chicken_one:
    """Map one library to reference genome using bowtie2

    Output SAM file is piped to samtools sort to generate a CRAM file.
    """
    input:
        forward_=BOWTIE2_NONHUMAN / "{sample}.{library}_1.fq.gz",
        reverse_=BOWTIE2_NONHUMAN / "{sample}.{library}_2.fq.gz",
        mock=BOWTIE2_HOSTS / "chicken_index",
        reference=REFERENCE / "chicken.fa.gz",
    output:
        cram=BOWTIE2_CHICKEN / "{sample}.{library}.cram",
    log:
        BOWTIE2_CHICKEN / "{sample}.{library}.log",
    benchmark:
        BOWTIE2_CHICKEN / "{sample}.{library}.bmk"
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
        runtime=6 * 60,
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


rule bowtie2_hosts_extract_nonchicken_one:
    """
    Keep only pairs unmapped to the human reference genome, sort by name rather
    than by coordinate, and convert to FASTQ.
    """
    input:
        cram=BOWTIE2_CHICKEN / "{sample}.{library}.cram",
        reference=REFERENCE / "chicken.fa.gz",
    output:
        forward_=BOWTIE2_NONCHICKEN / "{sample}.{library}_1.fq.gz",
        reverse_=BOWTIE2_NONCHICKEN / "{sample}.{library}_2.fq.gz",
    log:
        BOWTIE2_NONCHICKEN / "{sample}.{library}.log",
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
        | samtools sort \
            -n \
            -u \
            --threads {threads} \
        | samtools fastq \
            -1 {output.forward_} \
            -2 {output.reverse_} \
            -0 /dev/null \
            -c 9 \
            --threads {threads} \
        ) 2> {log} 1>&2
        """


rule bowtie2_hosts_extract_nonchicken_all:
    """Run bowtie2_extract_nonchicken_one for all libraries"""
    input:
        [
            BOWTIE2_NONCHICKEN / f"{sample}.{library}_{end}.fq.gz"
            for sample, library in SAMPLE_LIB
            for end in ["1", "2"]
        ],


rule bowtie2_hosts_report_all:
    """Generate bowtie2 reports for all libraries:
    - samtools stats
    - samtools flagstats
    - samtools idxstats
    """
    input:
        [
            f"{folder}/{sample}.{library}.{report}"
            for sample, library in SAMPLE_LIB
            for report in BAM_REPORTS
            for folder in [BOWTIE2_HUMAN, BOWTIE2_CHICKEN]
        ],


rule bowtie2_host:
    """Run bowtie2 on all libraries and generate reports"""
    input:
        rules.bowtie2_hosts_report_all.input,
        rules.bowtie2_hosts_extract_nonchicken_all.input,
