include: "functions.smk"
include: "nonpareil.smk"
include: "singlem.smk"
include: "coverm.smk"
include: "kraken2.smk"


rule stats:
    """Run the stats analyses: nonpareil and coverm"""
    input:
        rules.stats_nonpareil.output,
        # rules.stats_singlem.input,
        rules.stats_coverm.input,
        rules.stats_kraken2.input,


rule stats_with_singlem:
    """Run the nonpareil, coverm and singlem"""
    input:
        rules.stats.input,
        rules.stats_singlem.input,


localrules:
    stats_nonpareil,
    stats_singlem,
