#!/usr/bin/env Rscript
## Computing-time study: ICDT with M = 10 thinning repeats, the authors' default (10 x 8 + 8 fits)
## Run:    Rscript run_time_ICDT10.R                  # 6 (n,p) configurations x 100 replicates -> results/ICDT10/
##         TIME_WORKERS=10 Rscript run_time_ICDT10.R  # parallel single-threaded workers, one forked process per task
## Resume: run the same command again; finished tasks are skipped.
## Tables: Rscript summarize_time.R
## The script locates its own folder, so the whole time/ directory can be moved as is.
tc_script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else normalizePath(getwd())
}
ROOT <- tc_script_dir()
source(file.path(ROOT, "common", "timing_common.R"))
tc_main(list(root = ROOT, method = "ICDT10", family = "ICDT", K = NA_integer_, M = 10L))
