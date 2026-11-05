#!/usr/bin/env Rscript
## Computing-time study: TCV with K = 2 folds (TCV, TPCV, TICV, TIPCV share the 2 x 8 fits)
## Run:    Rscript run_time_TCV2.R                  # 6 (n,p) configurations x 100 replicates -> results/TCV2/
##         TIME_WORKERS=10 Rscript run_time_TCV2.R  # parallel single-threaded workers, one forked process per task
## Resume: run the same command again; finished tasks are skipped.
## Tables: Rscript summarize_time.R
## The script locates its own folder, so the whole time/ directory can be moved as is.
tc_script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else normalizePath(getwd())
}
ROOT <- tc_script_dir()
source(file.path(ROOT, "common", "timing_common.R"))
tc_main(list(root = ROOT, method = "TCV2", family = "TCV", K = 2L, M = NA_integer_))
