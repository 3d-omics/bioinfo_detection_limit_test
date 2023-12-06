# fastp ----
def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    adapter = samples[
        (samples["sample_id"] == wildcards.sample)
        & (samples["library_id"] == wildcards.library)
    ]["forward_adapter"].tolist()[0]
    if pd.isna(adapter):
        return "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    return adapter


def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    adapter = samples[
        (samples["sample_id"] == wildcards.sample)
        & (samples["library_id"] == wildcards.library)
    ]["reverse_adapter"].tolist()[0]
    if pd.isna(adapter):
        return "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
    return adapter


# def get_fastp_inputs(wildcards):
#     """Compose the inputs for fastp, depending on the library type"""
#     sample = wildcards.sample
#     library = wildcards.library
#     fastp_input_dict = {}
#     fastp_input_dict["forward_"] = READS / f"{sample}.{library}_1.fq.gz"
#     fastp_input_dict["reverse_"] = READS / f"{sample}.{library}_2.fq.gz"
#     return fastp_input_dict


# def get_fastp_outputs(wildcards):
#     """Compose the output for fastp, depending on the library type"""
#     sample = wildcards.sample
#     library = wildcards.library
#     fastp_output_dict = {}
#     fastp_output_dict["forward_"] = FASTP / f"{sample}.{library}_1.fq.gz"
#     fastp_output_dict["reverse_"] = FASTP / f"{sample}.{library}_2.fq.gz"
#     fastp_output_dict["unpaired1"] = FASTP / f"{sample}.{library}_u1.fq.gz"
#     fastp_output_dict["unpaired2"] = FASTP / f"{sample}.{library}_u2.fq.gz"
#     fastp_output_dict["html"] = FASTP / f"{sample}.{library}_fastp.html"
#     fastp_ouptut_dict["json"] = FASTP / f"{sample}.{library}_fastp.json"
#     return fastp_output_dict


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
        PRE_BOWTIE2
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
        PRE_BOWTIE2
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
    return PRE_BOWTIE2 / f"non{genome}" / f"{sample}.{library}_1.fq.gz"


def get_input_reverse_for_mag_mapping(wildcards):
    """Get the forward input file for mapping to MAGs"""
    sample = wildcards.sample
    library = wildcards.library
    if len(HOST_NAMES) == 0:
        return FASTP / f"{sample}.{library}_2.fq.gz"
    genome = HOST_NAMES[-1]
    return PRE_BOWTIE2 / f"non{genome}" / f"{sample}.{library}_2.fq.gz"


# def compose_input_string_for_bowtie2_mags_map_one_library_to_one_catalogue(wildcards):
#     forward_fn = get_input_forward_for_mag_mapping(wildcards)
#     reverse_fn = get_input_reverse_for_mag_mapping(wildcards)
#     if is_paired(wildcards):
#         return f"-1 {forward_fn} -2 {reverse_fn}"
#     return f"-U {forward_fn}"


# def compose_input_string_for_bowtie2_hosts_map_one(wildcards):
#     forward_fn = get_input_forward_for_host_mapping(wildcards)
#     reverse_fn = get_input_reverse_for_host_mapping(wildcards)
#     if is_paired(wildcards):
#         return f"-1 {forward_fn} -2 {reverse_fn}"
#     return f"-U {forward_fn}"


# def compose_rmdup_string_for_bowtie2_hosts_map_one(wildcards):
#     if is_paired(wildcards):
#         return ""
#     return "-s"


# def compose_output_string_for_bowtie2_hosts_extract_one(wildcards):
#     forward_fn = (
#         PRE_BOWTIE2
#         / f"non{wildcards.genome}"
#         / f"{wildcards.sample}.{wildcards.library}_1.fq.gz"
#     )
#     reverse_fn = (
#         PRE_BOWTIE2
#         / f"non{wildcards.genome}"
#         / f"{wildcards.sample}.{wildcards.library}_2.fq.gz"
#     )
#     if is_paired(wildcards):
#         return f"-1 {forward_fn} -2 {reverse_fn} -0 /dev/null"
#     return f"-0 {forward_fn}"


# def compose_filter_int_for_bowtie2_hosts_extract_one(wildcards):
#     return 12 if is_paired(wildcards) else 4


def get_kraken2_database(wildcards):
    """Get the kraken2 database"""
    return features["databases"]["kraken2"][wildcards.kraken2_db]


def compose_out_folder_for_pre_kraken2_assign_all(wildcards):
    """Compose the output folder for pre kraken2 assign all"""
    return KRAKEN2 / f"{wildcards.kraken2_db}"


def get_host_clean_forward(wildcards):
    """Get the forward input file that is clean from hosts"""
    last_genome = HOST_NAMES[-1]
    if len(HOST_NAMES) == 0:
        return FASTP / f"{wildcards.sample}.{wildcards.library}_1.fq.gz"
    return (
        PRE_BOWTIE2
        / f"non{last_genome}"
        / f"{wildcards.sample}.{wildcards.library}_1.fq.gz"
    )


def get_host_clean_reverse(wildcards):
    """Get the forward input file that is clean from hosts"""
    last_genome = HOST_NAMES[-1]
    if len(HOST_NAMES) == 0:
        return FASTP / f"{wildcards.sample}.{wildcards.library}_2.fq.gz"
    return (
        PRE_BOWTIE2
        / f"non{last_genome}"
        / f"{wildcards.sample}.{wildcards.library}_2.fq.gz"
    )


def compose_prefix_for_nonpareil(wildcards):
    """Compose prefix for nonpareil output files"""
    return NONPAREIL / f"{wildcards.sample}.{wildcards.library}"


def get_input_string_for_stats_singlem_pipe_one(wildcards):
    """Get the input string for stats singlem pipe one"""
    forward_fn = get_input_forward_for_stats(wildcards)
    reverse_fn = get_input_reverse_for_stats(wildcards)
    if is_paired(wildcards):
        return f"--forward {forward_fn} --reverse {reverse_fn}"
    else:
        return f"--forward {forward_fn}"
