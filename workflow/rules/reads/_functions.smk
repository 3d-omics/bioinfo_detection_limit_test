def get_forward(wildcards):
    """Get the forward read for a given sample and library"""
    forward_filename = samples.loc[
        (samples["sample_id"] == wildcards.sample)
        & (samples["library_id"] == wildcards.library)
    ]["forward_filename"].values[0]
    return forward_filename


def get_reverse(wildcards):
    """Get the reverse read for a given sample and library"""
    reverse_filename = samples.loc[
        (samples["sample_id"] == wildcards.sample)
        & (samples["library_id"] == wildcards.library)
    ]["reverse_filename"].values[0]
    return reverse_filename
