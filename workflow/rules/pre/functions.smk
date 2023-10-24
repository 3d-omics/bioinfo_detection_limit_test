# fastp ----
def get_forward_adapter(wildcards):
    """Get forward adapter for a sample and library."""
    return samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["forward_adapter"].tolist()[0]


def get_reverse_adapter(wildcards):
    """Get reverse adapter for a sample and library."""
    return samples[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["reverse_adapter"].tolist()[0]


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
        / f"non{prev_genome}/{wildcards.sample}.{wildcards.library}_1.fq.gz"
    )


def get_input_reverse_for_host_mapping(wildcards):
    """Get the reverse input file"""
    if wildcards.genome == HOST_NAMES[0]:
        return FASTP / f"{wildcards.sample}.{wildcards.library}_2.fq.gz"
    genome_index = HOST_NAMES.index(wildcards.genome)
    prev_genome = HOST_NAMES[genome_index - 1]
    return (
        BOWTIE2_HOSTS
        / f"non{prev_genome}/{wildcards.sample}.{wildcards.library}_2.fq.gz"
    )


def get_input_forward_for_mag_mapping(wildcards):
    """Get the forward input file for mapping to MAGs"""
    sample = wildcards.sample
    library = wildcards.library
    if len(HOST_NAMES) == 0:
        return FASTP / f"{sample}.{library}_1.fq.gz"
    genome = HOST_NAMES[-1]
    return BOWTIE2_HOSTS / f"non{genome}/{sample}.{library}_1.fq.gz"


def get_input_reverse_for_mag_mapping(wildcards):
    """Get the forward input file for mapping to MAGs"""
    sample = wildcards.sample
    library = wildcards.library
    if len(HOST_NAMES) == 0:
        return FASTP / f"{sample}.{library}_2.fq.gz"
    genome = HOST_NAMES[-1]
    return BOWTIE2_HOSTS / f"non{genome}/{sample}.{library}_2.fq.gz"
