RESULTS = Path("results/")


# All the other preprocessing folders are inherited from mg_preprocess
PRE_BOWTIE2 = RESULTS / "preprocess" / "bowtie2"

# folders for this pipeline
QUANT = Path("results/quantify/")
MAGS = QUANT / "mags"
QUANT_INDEX = QUANT / "build"
QUANT_BOWTIE2 = QUANT / "bowtie2"
COVERM = QUANT / "coverm"
