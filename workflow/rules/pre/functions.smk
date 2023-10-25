# fastp ----
def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    adapter = samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["forward_adapter"].tolist()[0]
    if pd.isna(adapter):
        return "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    return adapter


def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    adapter = samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["reverse_adapter"].tolist()[0]
    if pd.isna(adapter):
        return "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    return adapter


def get_fastp_inputs(wildcards):
    """Compose the inputs for fastp, depending on the library type"""
    sample = wildcards.sample
    library = wildcards.library
    sample_is_paired = is_paired(wildcards)
    fastp_input_dict = {}
    fastp_input_dict["forward_"] = READS / f"{sample}.{library}_1.fq.gz"
    fastp_input_dict["reverse_"] = (
        READS / f"{sample}.{library}_2.fq.gz" if sample_is_paired else "/dev/null"
    )
    return fastp_input_dict


def get_fastp_outputs(wildcards):
    """Compose the output for fastp, depending on the library type"""
    sample = wildcards.sample
    library = wildcards.library
    sample_is_paired = is_paired(wildcards)
    fastp_output_dict = {}
    fastp_output_dict["forward_"] = FASTP / f"{sample}.{library}_1.fq.gz"
    fastp_output_dict["reverse_"] = (
        FASTP / f"{sample}.{library}_2.fq.gz" if sample_is_paired else ""
    )
    fastp_output_dict["unpaired1"] = (
        FASTP / f"{sample}.{library}_u1.fq.gz" if sample_is_paired else ""
    )
    fastp_output_dict["unpaired2"] = (
        FASTP / f"{sample}.{library}_u2.fq.gz" if sample_is_paired else ""
    )
    fastp_output_dict["html"] = FASTP / f"{sample}.{library}_fastp.html"
    fastp_ouptut_dict["json"] = FASTP / f"{sample}.{library}_fastp.json"
    return fastp_output_dict


def compose_input_string_for_fastp_trim_one(wildcards):
    forward_fn = READS / f"{wildcards.sample}.{wildcards.library}_1.fq.gz"
    reverse_fn = READS / f"{wildcards.sample}.{wildcards.library}_2.fq.gz"
    if is_paired(wildcards):
        return f"--in1 {forward_fn} --in2 {reverse_fn}"
    else:
        return f"--in1 {forward_fn}"


def compose_output_string_for_fastp_trim_one(wildcards, threads):
    bgzip_string = lambda filename, threads: f">(bgzip -l 1 -@ {threads} > {filename})"
    forward_fn = bgzip_string(
        FASTP / f"{wildcards.sample}.{wildcards.library}_1.fq.gz", threads
    )
    reverse_fn = bgzip_string(
        FASTP / f"{wildcards.sample}.{wildcards.library}_2.fq.gz", threads
    )
    unpaired1_fn = bgzip_string(
        FASTP / f"{wildcards.sample}.{wildcards.library}_u1.fq.gz", threads
    )
    unpaired2_fn = bgzip_string(
        FASTP / f"{wildcards.sample}.{wildcards.library}_u2.fq.gz", threads
    )
    if is_paired(wildcards):
        return f"--out1 {forward_fn} --out2 {reverse_fn} --unpaired1 {unpaired1_fn} --unpaired2 {unpaired2_fn}"
    else:
        return f"--out1 {forward_fn}"


def compose_adapter_string_for_fastp_trim_one(wildcards):
    adapter_forward = get_forward_adapter(wildcards)
    adapter_reverse = get_reverse_adapter(wildcards)
    if is_paired(wildcards):
        return f"--adapter_sequence {adapter_forward} --adapter_sequence_r2 {adapter_reverse}"
    else:
        return f"--adapter_sequence {adapter_forward}"


# bowtie2 ----
def compose_rg_id(wildcards):
    """Compose the read group ID for bowtie2"""
    return f"{wildcards.sample}_{wildcards.library}"


def compose_rg_extra(wildcards):
    """Compose the read group extra information for bowtie2"""
    return f"LB:truseq_{wildcards.library}\tPL:Illumina\tSM:{wildcards.sample}"


def get_input_forward_for_host_mapping(wildcards):
    """Compose the forward input file"""
    if wildcards.genome == HOST_NAMES[0]:
        return FASTP / f"{wildcards.sample}.{wildcards.library}_1.fq.gz"
    genome_index = HOST_NAMES.index(wildcards.genome)
    prev_genome = HOST_NAMES[genome_index - 1]
    return (
        BOWTIE2_HOSTS
        / f"non{prev_genome}"
        / f"{wildcards.sample}.{wildcards.library}_1.fq.gz"
    )


def get_input_reverse_for_host_mapping(wildcards):
    """Get the reverse input file"""
    if wildcards.genome == HOST_NAMES[0]:
        return FASTP / f"{wildcards.sample}.{wildcards.library}_2.fq.gz"
    genome_index = HOST_NAMES.index(wildcards.genome)
    prev_genome = HOST_NAMES[genome_index - 1]
    return (
        BOWTIE2_HOSTS
        / f"non{prev_genome}"
        / f"{wildcards.sample}.{wildcards.library}_2.fq.gz"
    )


def get_input_forward_for_mag_mapping(wildcards):
    """Get the forward input file for mapping to MAGs"""
    sample = wildcards.sample
    library = wildcards.library
    if len(HOST_NAMES) == 0:
        return FASTP / f"{sample}.{library}_1.fq.gz"
    genome = HOST_NAMES[-1]
    return BOWTIE2_HOSTS / f"non{genome}" / f"{sample}.{library}_1.fq.gz"


def get_input_reverse_for_mag_mapping(wildcards):
    """Get the forward input file for mapping to MAGs"""
    sample = wildcards.sample
    library = wildcards.library
    if len(HOST_NAMES) == 0:
        return FASTP / f"{sample}.{library}_2.fq.gz"
    genome = HOST_NAMES[-1]
    return BOWTIE2_HOSTS / f"non{genome}" / f"{sample}.{library}_2.fq.gz"


def compose_input_string_for_bowtie2_mags_map_one_library_to_one_catalogue(wildcards):
    forward_fn = get_input_forward_for_mag_mapping(wildcards)
    reverse_fn = get_input_reverse_for_mag_mapping(wildcards)
    if is_paired(wildcards):
        return f"-1 {forward_fn} -2 {reverse_fn}"
    return f"-U {forward_fn}"


def compose_input_string_for_bowtie2_hosts_map_one(wildcards):
    forward_fn = get_input_forward_for_host_mapping(wildcards)
    reverse_fn = get_input_reverse_for_host_mapping(wildcards)
    if is_paired(wildcards):
        return f"-1 {forward_fn} -2 {reverse_fn}"
    return f"-U {forward_fn}"


def compose_rmdup_string_for_bowtie2_hosts_map_one(wildcards):
    if is_paired(wildcards):
        return ""
    return "-s"
