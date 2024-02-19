rule _preprocess__kraken2_database:
    input:
        database=get_kraken2_database,
    output:
        service(directory("/dev/shm/{kraken2_db}")),
    log:
        KRAKEN2 / "{kraken2_db}.log",
    conda:
        "__environment__.yml"
    resources:
        mem_mb=params["preprocess"]["kraken2"]["memory_gb"] * 1024,
    group:
        "preprocess__kraken2"
    shell:
        """
        mkdir --parents {output}

        rsync \
            --archive \
            --progress \
            --recursive \
            --times \
            --verbose \
            {input.database}/*.k2d \
            {output}/ \
        2> {log} 1>&2
        """


rule _preprocess__kraken2__assign:
    """Run kraken2 over one sample and using the database as a service"""
    input:
        forward_=FASTP / "{sample_id}.{library_id}_1.fq.gz",
        reverse_=FASTP / "{sample_id}.{library_id}_2.fq.gz",
        database="/dev/shm/{kraken2_db}",
    output:
        out_gz=KRAKEN2 / "{kraken2_db}" / "{sample_id}.{library_id}.out.gz",
        report=KRAKEN2 / "{kraken2_db}" / "{sample_id}.{library_id}.report",
    log:
        KRAKEN2 / "{kraken2_db}" / "{sample_id}.{library_id}.log",
    conda:
        "__environment__.yml"
    resources:
        mem_mb=1024,
        runtime=60,
    group:
        "preprocess__kraken2"
    shell:
        """
        kraken2 \
            --db {input.database} \
            --gzip-compressed \
            --paired \
            --output >(gzip > {output.out_gz}) \
            --report {output.report} \
            --memory-mapping \
            {input.forward_} \
            {input.reverse_} \
        2> {log} 1>&2
        """


rule preprocess__kraken2:
    """Run kraken2 over all samples at once using the /dev/shm/ trick."""
    input:
        [
            KRAKEN2 / kraken2_db / f"{sample_id}.{library_id}.report"
            for sample_id, library_id in SAMPLE_LIBRARY
            for kraken2_db in KRAKEN2_DBS
        ],
