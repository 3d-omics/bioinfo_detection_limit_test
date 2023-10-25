READS = Path("results/reads/")
REFERENCE = Path("results/reference/")

PRE = Path("results/preprocessing/")
FASTP = PRE / "fastp"
BOWTIE2_HOSTS = PRE / "bowtie2_host"
BOWTIE2_MAGS = PRE / "bowtie2_mags"


# BOWTIE2 = Path("results/bowtie2/")

STATS = Path("results/stats/")
KRAKEN2 = STATS / "kraken2"
NONPAREIL = STATS / "nonpareil"
SINGLEM = STATS / "singlem"
COVERM = STATS / "coverm"

REPORT_STEP = Path("reports/by_step/")
REPORT_LIBRARY = Path("reports/by_library/")
