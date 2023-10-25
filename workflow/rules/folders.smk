READS = Path("results/reads/")
REFERENCE = Path("results/reference/")

PRE = Path("results/preprocessing/")
FASTP = PRE / "fastp"
BOWTIE2_HOSTS = PRE / "bowtie2_host"
BOWTIE2_MAGS = PRE / "bowtie2_mags"

STATS = Path("results/stats/")
KRAKEN2 = STATS / "kraken2"
NONPAREIL = STATS / "nonpareil"
SINGLEM = STATS / "singlem"
COVERM = STATS / "coverm"

REPORT = Path("reports")
REPORT_STEP = REPORT / "by_step"
REPORT_LIBRARY = REPORT / "by_library"
