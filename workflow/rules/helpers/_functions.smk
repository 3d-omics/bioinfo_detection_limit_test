def get_sample_lib_pe():
    """Get all paired-end sample-libraries"""
    sample_lib_pe = samples[samples["library_type"] == "PE"][
        ["sample", "library"]
    ].values.tolist()
    return sample_lib_pe


def get_sample_lib_pe():
    """Get all single-end sample-libraries"""
    sample_lib_pe = get_sample_lib_pe()
    sample_lib_se = [
        samples[sample_library]
        for sample_library in sample_lib
        if sample_library not in sample_lib_pe
    ]
    return sample_lib_se


def is_paired(wildcards):
    """Test if a sample-library is paired-end"""
    sample_lib = [wildcards.sample, wildcards.library]
    if sample_lib in SAMPLE_LIB_PE:
        return True
    elif sample_lib in SAMPLE_LIB_SE:
        return False
    else:
        exit("Sample-library combination not found in config file.")


def double_ram(initial_memory):
    """Double the memory for each attempt"""
    return lambda wildcards, attempt: initial_memory * 2 ** (attempt - 1) * 1024
