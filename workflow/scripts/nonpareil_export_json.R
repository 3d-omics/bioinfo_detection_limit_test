#!/usr/bin/env Rscript

library(tidyverse)
library(argparse)

parser <- ArgumentParser()

parser$add_argument(
  "-i", "--input-npo",
  type = "character",
  dest = "input_npo",
  help = "path to Nonpareil's npo file"
)

parser$add_argument(
  "-j", "--output-json",
  type = "character",
  dest = "output_json",
  help = "Output JSON file"
)

# copied from https://snakemake-wrappers.readthedocs.io/en/stable/wrappers/nonpareil/plot.html#code
export_curve <- function(object) {
  # Extract variables
  n <- names(attributes(object))[c(1:12, 21:29)]
  x <- sapply(n, function(v) attr(object, v))
  names(x) <- n

  # Extract vectors
  n <- names(attributes(object))[13:20]
  y <- lapply(n, function(v) attr(object, v))
  names(y) <- n
  curve_json <- c(x, y)

  # Add model
  if (object$has.model) {
    # https://github.com/lmrodriguezr/nonpareil/blob/162f1697ab1a21128e1857dd87fa93011e30c1ba/utils/Nonpareil/R/Nonpareil.R#L330-L332
    x_min <- 1e3
    x_max <- signif(tail(attr(object, "x.adj"), n = 1) * 10, 1)
    x.model <- exp(seq(log(x_min), log(x_max), length.out = 1e3))
    y.model <- predict(object, lr = x.model)
    curve_json <- append(curve_json, list(x.model = x.model))
    curve_json <- append(curve_json, list(y.model = y.model))
  }

  curve_json
}

export_set <- function(object) {
  y <- lapply(object$np.curves, "export_curve")
  names(y) <- sapply(object$np.curves, function(n) n$label)
  jsonlite::prettify(jsonlite::toJSON(y, auto_unbox = TRUE))
}


args <- parser$parse_args()
input_npo <- args$input_npo
output_json <- args$output_json

dir.create(dirname(output_json), showWarnings = FALSE, recursive = TRUE)

if (file.info(input_npo)$size > 0) {
  input_npo %>%
    Nonpareil::Nonpareil.set(plot = FALSE) %>%
    export_set() %>%
    write(output_json)
} else {
  file.create(output_json)
}
