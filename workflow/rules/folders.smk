READS = Path("results/reads/")

REFERENCE = Path("results/reference/")

# All the other preprocessing folders are inherited from mg_preprocess
PRE_BOWTIE2 = Path("results/preprocess/bowtie2/")

QUANT = Path("results/quantify/")
QUANT_INDEX = QUANT / "index"
QUANT_BOWTIE2 = QUANT / "bowtie2"
COVERM = QUANT / "coverm"

REPORT = Path("reports")
REPORT_STEP = REPORT / "by_step"
REPORT_LIBRARY = REPORT / "by_library"
