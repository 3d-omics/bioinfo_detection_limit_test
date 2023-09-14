def get_forward(wildcards):
    """Get the forward read for a given sample and library"""
    forward_filename = samples.loc[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["forward"].values
    return forward_filename


def get_reverse(wildcards):
    """Get the reverse read for a given sample and library"""
    reverse_filename = samples.loc[
        (samples["sample"] == wildcards.sample)
        & (samples["library"] == wildcards.library)
    ]["reverse"].values
    return reverse_filename if not pd.isna(reverse_filename) else "/dev/null"
