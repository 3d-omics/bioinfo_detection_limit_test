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


def get_input_single_for_host_mapping(wildcards):
    """Get the reverse input file"""
    if wildcards.genome == HOST_NAMES[0]:
        return FASTP / f"{wildcards.sample}.{wildcards.library}_se.fq.gz"
    genome_index = HOST_NAMES.index(wildcards.genome)
    prev_genome = HOST_NAMES[genome_index - 1]
    return (
        BOWTIE2_HOSTS
        / f"non{prev_genome}/{wildcards.sample}.{wildcards.library}_se.fq.gz"
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


def get_input_single_for_mag_mapping(wildcards):
    """Get the forward input file for mapping to MAGs"""
    sample = wildcards.sample
    library = wildcards.library
    if len(HOST_NAMES) == 0:
        return FASTP / f"{sample}.{library}_se.fq.gz"
    genome = HOST_NAMES[-1]
    return BOWTIE2_HOSTS / f"non{genome}/{sample}.{library}_se.fq.gz"
