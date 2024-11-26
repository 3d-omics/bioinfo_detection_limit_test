RESULTS = Path("results/")


# All the other preprocessing folders are inherited from mg_preprocess
PRE_BOWTIE2 = RESULTS / "preprocess" / "bowtie2"
PRE_CLEAN = RESULTS / "preprocess" / "clean"

# folders for this pipeline
QUANT = Path("results/quantify/")
QUANT_MAGS = QUANT / "mags"
QUANT_BUILD = QUANT / "build"
QUANT_BOWTIE2 = QUANT / "bowtie2"
QUANT_COVERM = QUANT / "coverm"
