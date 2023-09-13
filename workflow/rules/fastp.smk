rule fastp_pe_trim_one:
    """Run fastp on one PE library"""
    input:
        forward_=READS / "{sample}.{library}_1.fq.gz",
        reverse_=READS / "{sample}.{library}_2.fq.gz",
    output:
        forward_=temp(FASTP / "{sample}.{library}_1.fq.gz"),
        reverse_=temp(FASTP / "{sample}.{library}_2.fq.gz"),
        unpaired1=temp(FASTP / "{sample}.{library}_u1.fq.gz"),
        unpaired2=temp(FASTP / "{sample}.{library}_u2.fq.gz"),
        html=FASTP / "{sample}.{library}_pe_fastp.html",
        json=FASTP / "{sample}.{library}_pe_fastp.json",
    log:
        FASTP / "{sample}.{library}_pe.log",
    params:
        adapter_forward=get_forward_adapter,
        adapter_reverse=get_reverse_adapter,
        extra=params["fastp"]["extra"],
        length_required=params["fastp"]["length_required"],
    threads: 16
    resources:
        mem_mb=4 * 1024,
        runtime=240,
    conda:
        "../envs/fastp.yml"
    shell:
        """
        fastp \
            --in1 {input.forward_} \
            --in2 {input.reverse_} \
            --out1 >(bgzip -l 1 -@ {threads} > {output.forward_}) \
            --out2 >(bgzip -l 1 -@ {threads} > {output.reverse_}) \
            --unpaired1 >(bgzip -l 1 -@ {threads} > {output.unpaired1}) \
            --unpaired2 >(bgzip -l 1 -@ {threads} > {output.unpaired2}) \
            --html {output.html} \
            --json {output.json} \
            --compression 1 \
            --verbose \
            --adapter_sequence {params.adapter_forward} \
            --adapter_sequence_r2 {params.adapter_reverse} \
            --length_required {params.length_required} \
            --thread {threads} \
            {params.extra} \
        2> {log} 1>&2
        """


rule fastp_se_trim_one:
    """Run fastp on one SE library"""
    input:
        forward_=READS / "{sample}.{library}_se.fq.gz",
    output:
        forward_=temp(FASTP / "{sample}.{library}_se.fq.gz"),
        html=FASTP / "{sample}.{library}_se_fastp.html",
        json=FASTP / "{sample}.{library}_se_fastp.json",
    log:
        FASTP / "{sample}.{library}_se.log",
    params:
        adapter_forward=get_forward_adapter,
        extra=params["fastp"]["extra"],
        length_required=params["fastp"]["length_required"],
    threads: 16
    resources:
        mem_mb=4 * 1024,
        runtime=240,
    conda:
        "../envs/fastp.yml"
    shell:
        """
        fastp \
            --in1 {input.forward_} \
            --out1 >(bgzip -l 1 -@ {threads} > {output.forward_}) \
            --html {output.html} \
            --json {output.json} \
            --compression 1 \
            --verbose \
            --adapter_sequence {params.adapter_forward} \
            --length_required {params.length_required} \
            --thread {threads} \
            {params.extra} \
        2> {log} 1>&2
        """


rule fastp_trim_all:
    """Run fastp over all libraries"""
    input:
        [
            FASTP / f"{sample}.{library}_{end}.fq.gz"
            for sample, library in SAMPLE_LIB_PE
            for end in "1 2 u1 u2".split(" ")
        ],
        [
            FASTP / f"{sample}.{library}_{end}.fq.gz"
            for sample, library in SAMPLE_LIB_SE
            for end in ["se"]
        ],


rule fastp_fastqc_all:
    """Run fastqc over all libraries"""
    input:
        [
            FASTP / f"{sample}.{library}_{end}_fastqc.{extension}"
            for sample, library in SAMPLE_LIB_PE
            for end in ["1", "2"]
            for extension in ["html", "zip"]
        ],
        [
            FASTP / f"{sample}.{library}_{end}_fastqc.{extension}"
            for sample, library in SAMPLE_LIB_SE
            for end in ["se"]
            for extension in ["html", "zip"]
        ],


rule fastp_report_all:
    """Collect fastp and fastqc reports"""
    input:
        [
            FASTP / f"{sample}.{library}_pe_fastp.json"
            for sample, library in SAMPLE_LIB_PE
        ],
        [
            FASTP / f"{sample}.{library}_se_fastp.json"
            for sample, library in SAMPLE_LIB_SE
        ],
        rules.fastp_fastqc_all.input,


rule fastp:
    """Run fastp and collect reports"""
    input:
        rules.fastp_trim_all.input,
        rules.fastp_report_all.input,
