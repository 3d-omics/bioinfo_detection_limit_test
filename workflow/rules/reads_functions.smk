def get_forward(wildcards):
    """Get the forward read for a given sample and library"""
    return samples.loc[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["forward"]


def get_reverse(wildcards):
    """Get the reverse read for a given sample and library"""
    return samples.loc[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["reverse"]
