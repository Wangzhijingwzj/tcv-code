#!/usr/bin/env Rscript
## Build the comparison tables from whatever the six run scripts have saved so far.
## Usage: Rscript summarize_time.R   -> results/summary/{all_rows.csv, summary_by_config.csv, REPORT.md}
tc_script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else normalizePath(getwd())
}
ROOT <- tc_script_dir()
source(file.path(ROOT, "common", "timing_common.R"))
tc_summarize(ROOT)
