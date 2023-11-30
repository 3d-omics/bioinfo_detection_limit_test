include: "_functions.smk"
include: "nonpareil.smk"
include: "singlem.smk"
include: "coverm.smk"
include: "kraken2.smk"


rule stats:
    """Run the stats analyses: nonpareil and coverm"""
    input:
        rules.stats__nonpareil.input,
        # rules.stats_singlem.input,
        rules.stats__coverm.input,
        rules.stats__kraken2.input,


rule stats_with_singlem:
    """Run the nonpareil, coverm and singlem"""
    input:
        rules.stats.input,
        rules.stats__singlem.input,
