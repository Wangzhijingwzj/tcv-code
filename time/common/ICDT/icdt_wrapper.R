################################################################################
## icdt_wrapper.R
## Thin wrapper around the authors' ICDT code
##   Wang, Xu, Zhao & Wang (2025, JASA), "Data thinning for Poisson factor
##   models and its applications".  Upstream: https://github.com/Wangzhijingwzj/ICDT
##
## The authors' computational code is retained; only the sourceCpp path is adapted.
## This file only
##   (1) loads the packages that the upstream real-data script real_bbc.R loads,
##   (2) sources functions.R from common/ICDT under the runner's ROOT
##       (its sourceCpp path points to common/ICDT/bgfm.cpp),
##   (3) defines run_icdt(), which calls the upstream fac_ic() on one count
##       matrix and turns its criterion values into selected factor numbers
##       exactly as the upstream scripts do:
##         ICDT    = which.min(ICDT_val)          (ICDT.R line 144, real_bbc.R line 46)
##         IC4/IC5 = which.min(IC4_val / IC5_val) (ICDT.R lines 142-143)
##         VSDT    = floor(getModeAverage(which.min of each thinning))
##                                                (ICDT.R line 166, column "fis_B10")
##
## Usage
##   ROOT <- path.expand("~/multiDT/TCV-ICDT"); source(file.path(ROOT, "common/ICDT/icdt_wrapper.R"))
##   res <- run_icdt(x, qmax = 8, reptimes = 10, maxit = 30)
##   res$ICDT; res$VSDT; res$VSDT_each; res$IC4; res$IC5; res$time_sec
################################################################################

.icdt_need <- function(path, what) {
  if (!file.exists(path))
    stop(sprintf("%s not found at %s. Check the common/ICDT files under ROOT.",
                 what, path), call. = FALSE)
}
.icdt_need(file.path(ROOT, "common/ICDT/functions.R"), "upstream functions.R")
.icdt_need(file.path(ROOT, "common/ICDT/bgfm.cpp"), "upstream bgfm.cpp")

## Packages attached by the upstream real-data script (real_bbc.R).
## functions.R itself attaches Rcpp and RcppArmadillo and compiles bgfm.cpp (~5 s).
suppressPackageStartupMessages({
  library(irlba)
  library(doSNOW)
  library(MASS)
})
source(file.path(ROOT, "common/ICDT/functions.R"))   # computational code unchanged

## run_icdt: ICDT / VSDT / IC4 / IC5 for one n x p count matrix x
##   qmax     : candidate numbers of factors 1..qmax   (upstream simulations: 8)
##   reptimes : number of thinning repetitions M       (upstream: 10)
##   maxit    : maxIter passed to GFM::gfm             (upstream simulations: default 30)
##   quiet    : capture the console output of fac_ic() instead of printing it
##   restore_rng : fac_ic() calls set.seed(b + 1000) for every thinning b, which
##                 overwrites the caller's RNG stream; TRUE saves and restores
##                 .Random.seed so a simulation loop around run_icdt() is unaffected
run_icdt <- function(x, qmax = 8, reptimes = 10, maxit = 30,
                     quiet = TRUE, restore_rng = TRUE) {
  x <- as.matrix(x)
  storage.mode(x) <- "double"
  n <- nrow(x); p <- ncol(x)
  if (qmax < 2) stop("qmax must be at least 2")
  if (restore_rng) {
    had_seed <- exists(".Random.seed", envir = globalenv(), inherits = FALSE)
    if (had_seed) old_seed <- get(".Random.seed", envir = globalenv())
    on.exit({
      if (had_seed) {
        assign(".Random.seed", old_seed, envir = globalenv())
      } else if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) {
        rm(".Random.seed", envir = globalenv())
      }
    }, add = TRUE)
  }
  t0 <- proc.time()[["elapsed"]]
  if (quiet) {
    console <- utils::capture.output(
      re <- fac_ic(x, qmax = qmax, reptimes = reptimes, maxit = maxit),
      type = "output")
  } else {
    re <- fac_ic(x, qmax = qmax, reptimes = reptimes, maxit = maxit)
    console <- character(0)
  }
  elapsed <- proc.time()[["elapsed"]] - t0
  vsdt_each <- apply(re$VSDT_val_mat, 1, which.min)       # one VSDT estimate per thinning
  f_np <- log(n * p / (n + p))                             # penalty f(n,p) used upstream
  cov_vec <- (re$ICDT_val - re$loglike_val) / f_np         # = colMeans(cov_mat) inside fac_ic
  list(
    ICDT        = which.min(re$ICDT_val),
    VSDT        = floor(getModeAverage(vsdt_each)),        # "fis_B10" of ICDT.R
    VSDT_each   = vsdt_each,
    IC4         = which.min(re$IC4_val),
    IC5         = which.min(re$IC5_val),
    cov_penalty = cov_vec,                                 # per-entry scale, as inside fac_ic
    df_hat      = cov_vec * n * p,                         # count scale; compare with (r+1)p + rn
    values      = re,                                      # raw output of upstream fac_ic()
    n = n, p = p, qmax = qmax, reptimes = reptimes, maxit = maxit,
    n_fits      = reptimes * qmax + qmax,                  # GFM fits performed
    time_sec    = unname(elapsed),
    console     = console
  )
}
