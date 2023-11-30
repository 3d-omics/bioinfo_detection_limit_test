rule _preprocess__fastp__trim:
    """Run fastp on one PE library"""
    input:
        unpack(get_fastp_inputs),
    output:
        forward_=temp(FASTP / "{sample}.{library}_1.fq.gz"),
        reverse_=temp(touch(FASTP / "{sample}.{library}_2.fq.gz")),
        unpaired1=temp(touch(FASTP / "{sample}.{library}_u1.fq.gz")),
        unpaired2=temp(touch(FASTP / "{sample}.{library}_u2.fq.gz")),
        html=FASTP / "{sample}.{library}_fastp.html",
        json=FASTP / "{sample}.{library}_fastp.json",
    log:
        FASTP / "{sample}.{library}.log",
    params:
        extra=params["pre"]["fastp"]["extra"],
        length_required=params["pre"]["fastp"]["length_required"],
        input_string=compose_input_string_for_fastp_trim_one,
        output_string=compose_output_string_for_fastp_trim_one,
        adapter_string=compose_adapter_string_for_fastp_trim_one,
    threads: 16
    resources:
        mem_mb=4 * 1024,
        runtime=240,
    conda:
        "_env.yml"
    shell:
        """
        fastp \
            {params.input_string} \
            {params.output_string} \
            {params.adapter_string} \
            --html {output.html} \
            --json {output.json} \
            --compression 1 \
            --verbose \
            --length_required {params.length_required} \
            --thread {threads} \
            {params.extra} \
        2> {log} 1>&2
        """


rule preprocess__fastp__trim:
    """Run fastp over all libraries"""
    input:
        [
            FASTP / f"{sample}.{library}_{end}.fq.gz"
            for sample, library in SAMPLE_LIB
            for end in "1 2 u1 u2".split(" ")
        ],


rule preprocess__fastp__fastqc:
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
            for end in ["1"]
            for extension in ["html", "zip"]
        ],


rule preprocess__fastp__report:
    """Collect fastp and fastqc reports"""
    input:
        [FASTP / f"{sample}.{library}_fastp.json" for sample, library in SAMPLE_LIB],
        rules.preprocess__fastp__fastqc.input,


rule preprocess__fastp:
    """Run fastp and collect reports"""
    input:
        rules.preprocess__fastp__trim.input,
        rules.preprocess__fastp__report.input,
