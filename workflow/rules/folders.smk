READS = Path("results/reads/")
REFERENCE = Path("results/reference/")
FASTP = Path("results/fastp/")
KRAKEN2 = Path("results/kraken2/")

# BOWTIE2 = Path("results/bowtie2/")
BOWTIE2_HOSTS = Path("results/bowtie2_host")
BOWTIE2_HUMAN = BOWTIE2_HOSTS / "human"
BOWTIE2_NONHUMAN = BOWTIE2_HOSTS / "nonhuman"
BOWTIE2_CHICKEN = BOWTIE2_HOSTS / "chicken"
BOWTIE2_NONCHICKEN = BOWTIE2_HOSTS / "nonchicken"

BOWTIE2_MAGS = Path("results/bowtie2_mags")

STATS = Path("results/stats/")
STATS_NONPAREIL = STATS / "nonpareil"
STATS_SINGLEM = STATS / "singlem"
STATS_COVERM = STATS / "coverm"

REPORT_STEP = Path("reports/by_step/")
REPORT_LIBRARY = Path("reports/by_library/")
