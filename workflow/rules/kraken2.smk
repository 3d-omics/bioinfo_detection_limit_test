rule kraken2_assign_one:
    """Run kraken2 over one library and using one database."""
    input:
        forward_=FASTP / "{sample}.{library}_1.fq.gz",
        reverse_=FASTP / "{sample}.{library}_2.fq.gz",
        database=get_kraken2_database,
    output:
        out_gz=KRAKEN2 / "{kraken2_db}/{sample}.{library}.out.gz",
        report=KRAKEN2 / "{kraken2_db}/{sample}.{library}.report",
    log:
        log=KRAKEN2 / "{kraken2_db}/{sample}.{library}.log",
    conda:
        "../envs/kraken2.yml"
    threads: 24
    resources:
        mem_mb=params["kraken2"]["mem_mb"],
        runtime=60,
    shell:
        """
        kraken2 \
            --db {input.database} \
            --threads {threads} \
            --paired \
            --gzip-compressed \
            --output >(pigz > {output.out_gz}) \
            --report {output.report} \
            {input.forward_} \
            {input.reverse_} \
        > {log} 2>&1
        """


rule kraken2_assign_all:
    input:
        [
            KRAKEN2 / f"{kraken2_db}/{sample}.{library}.report"
            for sample, library in SAMPLE_LIB
            for kraken2_db in KRAKEN2_DBS
        ],


rule kraken2_report_one:
    input:
        KRAKEN2 / "{kraken2_db}/{sample}.{library}.report",


rule kraken2_report_all:
    input:
        [
            KRAKEN2 / f"{kraken2_db}/{sample}.{library}.report"
            for sample, library in SAMPLE_LIB
            for kraken2_db in KRAKEN2_DBS
        ],


rule kraken2:
    input:
        rules.kraken2_assign_all.input,
        rules.kraken2_report_all.input,
