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
