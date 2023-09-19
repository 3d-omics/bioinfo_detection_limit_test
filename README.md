# Snakemake workflow: `detection limit test`

[![Snakemake](https://img.shields.io/badge/snakemake-≥6.3.0-brightgreen.svg)](https://snakemake.github.io)
[![GitHub actions status](https://github.com/3d-omics/bioinfo_detection_limit_test/workflows/Tests/badge.svg?branch=devel)](https://github.com/3d-omics/bioinfo_detection_limit_test/actions?query=branch%devel+workflow%3ATests)


A Snakemake workflow for assessing detection limit from laser-microdissected samples.

## Usage

0. Requirements
   1.  [`miniconda`](https://docs.conda.io/en/latest/miniconda.html) / [`mamba`](https://mamba.readthedocs.io)
   2.  [`snakemake`](snakemake.readthedocs.io/)

1. Clone the repository
Clone the repository, and set it as the working directory.

```
git clone --recursive https://github.com/3d-omics/detection_limit_test.git
cd detection_limit_test
```

2. Run the pipeline with the test data (takes 5 minutes to download the required software)
```
snakemake \
    --use-conda \
    --conda-frontend mamba \
    -j 8
```

3. Edit the following files:
   1. `config/samples.tsv`: the control file with the sequencing libraries and their location.
   2. `config/features.tsv`: the references against which to map the libraries: human, chicken / pig, MAG catalogue.
   3. `config/params.tsv`: parameters for every program. The defaults are reasonable.

4. Run the pipeline and go for a walk:

```
./run  # locally
./run_slurm  # in a cluster environment with slurm
```

## Brief description

1. Trim reads and remove adaptors with `fastp`
2. Map to human, chicken / pig, mag catalogue:
   1. Map to the reference with `bowtie2`
   2. Extract the reads that have one of both ends unmapped with `samtools`
   3. Map those unmapped reads to the next reference
3. Generate MAG statistics with
   1. `coverm`
   2. `singlem`
   3. `nonpareil`
4. Generate lots of reports in the `reports/` folder

## Rulegraph

![rulegraph](rulegraph.svg)


## Possible problems

- `singlem` and/or `nonpareil` didnot finish some output because of low coverage.
  Paste this:

   ```
   Rscript workflow/scripts/aggregate_nonpareil.R \
      --input-folder results/stats/nonpareil \
      --output-file results/stats/nonpareil.tsv

   Rscript workflow/scripts/aggregate_singlem.R \
      --input-folder results/stats/singlem \
      --output-file results/stats/singlem.tsv
   ```

- Or comment with a `#` the `samples.tsv` file.


## References

- [fastp](https://github.com/OpenGene/fastp)
- [bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)
- [samtools](https://www.htslib.org/)
- [coverm](https://github.com/wwood/CoverM)
- [singlem](https://github.com/wwood/singlem)
- [nonpareil](http://enve-omics.ce.gatech.edu/nonpareil/)
- [fastqc](https://github.com/s-andrews/FastQC)
- [multiqc](https://multiqc.info/)
