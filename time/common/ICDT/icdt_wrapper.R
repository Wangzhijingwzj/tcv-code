################################################################################
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
