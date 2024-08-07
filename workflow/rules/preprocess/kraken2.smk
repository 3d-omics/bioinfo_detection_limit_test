rule preprocess__kraken2__assign__:
    """
    Run kraken2 over all samples at once using the /dev/shm/ trick.

    NOTE: /dev/shm may be not empty after the job is done.
    """
    input:
        forwards=[
            FASTP / f"{sample}.{library}_1.fq.gz" for sample, library in SAMPLE_LIBRARY
        ],
        rerverses=[
            FASTP / f"{sample}.{library}_2.fq.gz" for sample, library in SAMPLE_LIBRARY
        ],
        database=get_kraken2_database,
    output:
        out_gzs=[
            KRAKEN2 / "{kraken2_db}" / f"{sample}.{library}.out.gz"
            for sample, library in SAMPLE_LIBRARY
        ],
        reports=[
            KRAKEN2 / "{kraken2_db}" / f"{sample}.{library}.report"
            for sample, library in SAMPLE_LIBRARY
        ],
    log:
        KRAKEN2 / "{kraken2_db}.log",
    params:
        in_folder=FASTP,
        out_folder=compose_out_folder_for_eval_kraken2_assign_all,
        kraken_db_name="{kraken2_db}",
        samples=lambda w: [f"{sample}.{library}" for sample, library in SAMPLE_LIBRARY],
    conda:
        "__environment__.yml"
    shell:
        """
        {{
            echo Running kraken2 in $(hostname) 2>> {log} 1>&2

            mkdir --parents /dev/shm/{params.kraken_db_name}
            mkdir --parents {params.out_folder}

            rsync \
                --archive \
                --progress \
                --recursive \
                --times \
                --verbose \
                --chown $(whoami):$(whoami) \
                --chmod u+rw \
                {input.database}/*.k2d \
                /dev/shm/{params.kraken_db_name} \
            2>> {log} 1>&2

            parallel \
                --line-buffer \
                kraken2 \
                    --db /dev/shm/{params.kraken_db_name} \
                    --threads 1 \
                    --gzip-compressed \
                    --paired \
                    --output ">("gzip ">" {params.out_folder}/{{}}.out.gz")" \
                    --report {params.out_folder}/{{}}.report \
                    {params.in_folder}/{{}}_1.fq.gz \
                    {params.in_folder}/{{}}_2.fq.gz \
                "2>" {params.out_folder}/{{}}.log  "1>&2" \
            ::: {params.samples}

        }} || {{

            echo "Failed job" 2>> {log} 1>&2

        }}

        rm \
            --force \
            --recursive \
            --verbose \
            /dev/shm/{params.kraken_db_name} \
        2>>{log} 1>&2
        """


rule preprocess__kraken2:
    """Run kraken2 over all samples at once using the /dev/shm/ trick."""
    input:
        [
            KRAKEN2 / kraken2_db / f"{sample_id}.{library_id}.{extension}"
            for sample_id, library_id in SAMPLE_LIBRARY
            for kraken2_db in KRAKEN2_DBS
            for extension in ["out.gz", "report"]
        ],
