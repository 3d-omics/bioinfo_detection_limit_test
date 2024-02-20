READS = Path("results/reads/")

REFERENCE = Path("results/reference/")

PRE = Path("results/preprocess/")
FASTP = PRE / "fastp"
PRE_INDEX = PRE / "index"
PRE_BOWTIE2 = PRE / "bowtie2"
KRAKEN2 = PRE / "kraken2"
NONPAREIL = PRE / "nonpareil"
SINGLEM = PRE / "singlem"

QUANT = Path("results/quantify/")
QUANT_INDEX = QUANT / "index"
QUANT_BOWTIE2 = QUANT / "bowtie2"
COVERM = QUANT / "coverm"

REPORT = Path("reports")
REPORT_STEP = REPORT / "by_step"
REPORT_LIBRARY = REPORT / "by_library"
