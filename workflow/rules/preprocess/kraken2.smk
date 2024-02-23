rule preprocess__kraken2__:
    """
    Run kraken2 over all samples at once using the /dev/shm/ trick.

    NOTE: /dev/shm may be not empty after the job is done.
    """
    input:
        forwards=[
            FASTP / f"{sample_id}.{library_id}_1.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        rerverses=[
            FASTP / f"{sample_id}.{library_id}_2.fq.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        database=lambda w: features["databases"]["kraken2"][w.kraken_db],
    output:
        out_gzs=[
            KRAKEN2 / "{kraken_db}" / f"{sample_id}.{library_id}.out.gz"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
        reports=[
            KRAKEN2 / "{kraken_db}" / f"{sample_id}.{library_id}.report"
            for sample_id, library_id in SAMPLE_LIBRARY
        ],
    log:
        KRAKEN2 / "{kraken_db}.log",
    params:
        in_folder=FASTP,
        out_folder=lambda w: KRAKEN2 / w.kraken_db,
        kraken_db_shm="/dev/shm/{kraken_db}",
    conda:
        "__environment__.yml"
    shell:
        """
        {{
            echo Running kraken2 in $(hostname) 2>> {log} 1>&2

            mkdir --parents {params.kraken_db_shm}
            mkdir --parents {params.out_folder}

            rsync \
                --archive \
                --progress \
                --recursive \
                --times \
                --verbose \
                --chown $USER:$USER \
                --chmod u+rw \
                {input.database}/*.k2d \
                {params.kraken_db_shm} \
            2>> {log} 1>&2


            parallel \
                --jobs {threads} \
                kraken2 \
                    --db {params.kraken_db_shm} \
                    --threads 1 \
                    --gzip-compressed \
                    --paired \
                    --output {params.out_folder}/{{}}.out.gz \
                    --report {params.out_folder}/{{}}.report \
                    --memory-mapping \
                    {params.in_folder}/{{}}_1.fq.gz \
                    {params.in_folder}/{{}}_2.fq.gz \
                "2> {params.out_folder}/{{}}.log 1>&2" \
            ::: $(ls -1 {input.forwards} | xargs -I "{{}}" basename {{}} _1.fq.gz) \
            2>> {log} 1>&2

        }} || {{
            echo "Failed job" 2>> {log} 1>&2
        }}

        rm --force --recursive --verbose {params.kraken_db_shm} 2>>{log} 1>&2
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
