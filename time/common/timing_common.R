################################################################################
## timing_common.R
## Shared code for the TCV vs ICDT computing-time study (Example 1 data).
## Sourced by the six run_time_<METHOD>.R scripts; not meant to be run directly.
##
## Nothing in the timed code is modified. This file only loads existing code and
## records the time of one call per data set:
##   * TCV:  common/250113function_tcv.R and
##           common/confirm_jmle_omp_poisson_with_intercept_missing.cpp are the
##           files used by the paper's simulations (code/simu/simu-*). They are
##           loaded whole, exactly as the simulation workers did:
##             source("250113function_tcv.R"); Rcpp::sourceCpp("confirm_jmle_...cpp")
##           Timed call: multiDT(x, K = K, rmax = 8).
##   * ICDT: common/ICDT/ is the copy already used by the TCV-ICDT server package
##           (authors' functions.R with its one sourceCpp path adapted, plus
##           icdt_wrapper.R). Loaded as in TCV-ICDT/<scenario>/run_icdt_<scenario>.R.
##           Timed call: run_icdt(x, qmax = 8, reptimes = M, maxit = 30).
##   * Data: the Example 1 generator block of common/source/250410_TCV_norm.R is
##           read and evaluated with the same code as TCV-ICDT/norm/run_icdt_norm.R.
##
## Time = the whole method call, from the count matrix going in to the criterion
## values coming out: thinning, every model fit (including any fallback the method
## itself performs), predicted intensities and held-out losses. Data generation is
## timed separately and excluded. The nominal number of fits (K*8 for TCV_K,
## (M+1)*8 for ICDT_M) is recorded as structural information only.
##
## Each run script writes only to results/<METHOD>/. The cross-method tables are
## built afterwards by summarize_time.R, so concurrent jobs never share a file.
################################################################################

## Single-threaded BLAS / OpenMP for every process started from here. The
## Slurm template exports the same variables before R starts.
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")

RMAX    <- 8L      # candidate numbers of factors 1..RMAX
MAXITER <- 30L     # maxit passed to run_icdt(); multiDT() hard-codes maxIter = 30 itself
TRUE_Q  <- 5L      # Example 1 true number of factors (recorded only)
DATA_SEED_OFFSET <- 1000L   # data seed = 1000 + s, as in the original simulations

## Design: 6 (n, p) configurations x S replicates per method. The order of
## SIZES is also the run order: balanced configurations from small to large,
## the two unbalanced ones last.
SIZES <- data.frame(n = c(100L, 200L, 500L, 1000L, 100L, 1000L),
                    p = c(100L, 200L, 500L, 1000L, 1000L, 100L))
S <- 100L

METHOD_ORDER <- c("TCV2", "TCV5", "TCV10", "ICDT3", "ICDT5", "ICDT10")

## Number of concurrent single-threaded worker processes (forked, one per task).
WORKERS <- suppressWarnings(as.integer(Sys.getenv("TIME_WORKERS", "1")))
if (is.na(WORKERS) || WORKERS < 1L) WORKERS <- 1L

## ---------------------------------------------------------------------------
## Paths and logging
## ---------------------------------------------------------------------------
tc_files <- function(root) list(
  core      = file.path(root, "common", "250113function_tcv.R"),
  cpp       = file.path(root, "common", "confirm_jmle_omp_poisson_with_intercept_missing.cpp"),
  wrapper   = file.path(root, "common", "ICDT", "icdt_wrapper.R"),
  functions = file.path(root, "common", "ICDT", "functions.R"),
  generator = file.path(root, "common", "source", "250410_TCV_norm.R"))

tc_paths <- function(root, method) {
  res <- file.path(root, "results")
  list(root = root, results = res,
       method_dir = file.path(res, method),
       rows = file.path(res, method, "rows"),
       log = file.path(res, method, "run.log"))
}

tc_log <- function(paths, msg) {
  line <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), msg)
  cat(line, "\n")
  cat(line, "\n", file = paths$log, append = TRUE)
}

## ---------------------------------------------------------------------------
## Dependencies: every library() call in the files that will be sourced, plus
## the packages used by namespace.
## ---------------------------------------------------------------------------
tc_libs_in <- function(files) {
  l <- unlist(lapply(files, readLines, warn = FALSE))
  l <- sub("#.*$", "", l)
  m <- regmatches(l, regexpr("library\\(([A-Za-z0-9.]+)\\)", l))
  unique(sub("library\\(([A-Za-z0-9.]+)\\)", "\\1", m))
}

tc_check_packages <- function(files) {
  need <- unique(c(tc_libs_in(files), "GFM", "Rcpp", "RcppArmadillo",
                   "irlba", "pracma", "parallel"))
  ok <- vapply(need, requireNamespace, logical(1), quietly = TRUE)
  if (!all(ok))
    stop("Missing R packages: ", paste(need[!ok], collapse = ", "), call. = FALSE)
  invisible(need)
}

## ---------------------------------------------------------------------------
## Load the existing code, unchanged
## ---------------------------------------------------------------------------
tc_load_tcv <- function(f) {
  ## Same two statements as the workers of the original simulation scripts.
  for (x in c(f$core, f$cpp)) if (!file.exists(x)) stop("Missing ", x, call. = FALSE)
  suppressPackageStartupMessages(source(f$core))   # whole file, into the global environment
  Rcpp::sourceCpp(f$cpp)                           # defines confirm_CJMLE_poisson_cpp()
  invisible(TRUE)
}

tc_load_icdt <- function(f, root) {
  ## Same as TCV-ICDT/<scenario>/run_icdt_<scenario>.R: ROOT, irlba, Rfast, wrapper.
  for (x in c(f$wrapper, f$functions)) if (!file.exists(x)) stop("Missing ", x, call. = FALSE)
  assign("ROOT", root, envir = globalenv())        # read by icdt_wrapper.R and functions.R
  suppressPackageStartupMessages({ library(irlba); library(Rfast) })
  source(f$wrapper)                                # sources the authors' functions.R
  invisible(TRUE)
}

## Copied from TCV-ICDT/norm/run_icdt_norm.R: read only the original generator
## block, without running its simulation loop. Three separate add_identifiability
## calls, their RNG consumption and the Poisson loops are kept intact.
tc_load_generator <- function(f) {
  for (x in c(f$generator, f$core)) if (!file.exists(x)) stop("Missing ", x, call. = FALSE)
  source_lines <- readLines(f$generator, warn = FALSE)
  first <- grep("^[[:space:]]*set.seed\\(s \\+ 1000\\)", source_lines)
  last <- grep("^[[:space:]]*all_results\\$x <- x", source_lines)
  stopifnot(length(first) == 1L, length(last) == 1L, first < last)
  GENERATION_EXPR <- parse(text = c("all_results <- list()", source_lines[first:last]))
  helper <- Filter(function(e) is.call(e) &&
                     as.character(e[[1]]) %in% c("<-", "=") &&
                     identical(e[[2]], as.name("add_identifiability")),
                   as.list(parse(f$core)))
  GEN_ENV <- new.env(parent = globalenv())
  eval(helper[[1]], envir = GEN_ENV)
  GEN_ENV$Diag <- pracma::Diag
  GEN_ENV$irlba <- irlba::irlba
  GEN_ENV$has_full_zero_row_or_col <- function(mat)
    any(rowSums(mat == 0) == ncol(mat)) || any(colSums(mat == 0) == nrow(mat))
  assign("GENERATION_EXPR", GENERATION_EXPR, envir = globalenv())
  assign("GEN_ENV", GEN_ENV, envir = globalenv())
  invisible(TRUE)
}

generate_data <- function(n, p, s) {
  e <- new.env(parent = GEN_ENV)
  e$n <- n; e$p <- p; e$s <- s
  eval(GENERATION_EXPR, envir = e)
  list(x = e$x, q = e$q)
}

## ---------------------------------------------------------------------------
## Result files
## ---------------------------------------------------------------------------
tc_task_path <- function(paths, n, p, s)
  file.path(paths$rows, sprintf("n%d_p%d_s%03d.rds", n, p, s))

tc_n_fits_nominal <- function(cfg)
  if (cfg$family == "TCV") cfg$K * RMAX else (cfg$M + 1L) * RMAX

## Write a result atomically: a reader never sees a half-written .rds file.
tc_save_row <- function(row, path) {
  tmp <- paste0(path, ".tmp")
  saveRDS(row, tmp)
  if (!file.rename(tmp, path)) stop("Could not rename ", tmp, " to ", path, call. = FALSE)
  invisible(path)
}

## Validate an existing result file against the current settings before it is
## accepted as finished. Unreadable or mismatching files stop the run; they are
## never skipped silently and never overwritten. Returns TRUE when the file was
## produced with the same number of workers as the current run.
tc_check_existing <- function(path, cfg, n, p, s) {
  old <- tryCatch(readRDS(path), error = function(e) NULL)
  if (is.null(old) || !is.list(old))
    stop("Unreadable result file ", path,
         " (interrupted write?). Move it away and rerun.", call. = FALSE)
  same <- identical(old$method, cfg$method) && identical(old$family, cfg$family) &&
    isTRUE(old$n == n) && isTRUE(old$p == p) && isTRUE(old$s == s) &&
    isTRUE(old$seed == s + DATA_SEED_OFFSET) &&
    identical(as.integer(old$K), as.integer(cfg$K)) &&
    identical(as.integer(old$M), as.integer(cfg$M)) &&
    isTRUE(old$rmax == RMAX) && isTRUE(old$maxiter == MAXITER)
  if (!same)
    stop("Result file ", path, " does not match the current settings ",
         "(method, K/M, seed, rmax or maxIter). Move it away and rerun.", call. = FALSE)
  isTRUE(old$workers == WORKERS)
}

## Interim statistics: once every replicate of an (n, p) configuration has a
## result file, log mean / median / sd of the ok times of that configuration.
## Computed from the saved files after the task is finished, so the timing
## itself is untouched. (Two workers finishing the same configuration at the
## same moment may both log it; the lines are then identical.)
tc_config_progress <- function(paths, cfg, n, p) {
  files <- vapply(seq_len(S), function(s) tc_task_path(paths, n, p, s), "")
  if (!all(file.exists(files))) return(invisible(FALSE))
  d <- tc_read_rows(files)
  if (is.null(d)) return(invisible(FALSE))
  ok <- d[d$status == "ok", ]
  stats <- if (nrow(ok)) sprintf("mean=%.1f median=%.1f sd=%.1f min=%.1f max=%.1f seconds",
                                 mean(ok$elapsed_sec), median(ok$elapsed_sec),
                                 if (nrow(ok) > 1) sd(ok$elapsed_sec) else NA_real_,
                                 min(ok$elapsed_sec), max(ok$elapsed_sec)) else "no ok task"
  tc_log(paths, sprintf("CONFIG %-6s n=%4d p=%4d finished: ok=%d/%d %s",
                        cfg$method, n, p, nrow(ok), S, stats))
  invisible(TRUE)
}

## ---------------------------------------------------------------------------
## One task = one (n, p, s): generate data, time the method call, save a row
## ---------------------------------------------------------------------------
tc_run_task <- function(task, cfg, paths) {
  n <- as.integer(task$n); p <- as.integer(task$p); s <- as.integer(task$s)
  path <- tc_task_path(paths, n, p, s)
  if (file.exists(path)) return(invisible(path))   # validated in tc_main()
  started <- Sys.time()
  tc_log(paths, sprintf("START %-6s n=%4d p=%4d s=%3d", cfg$method, n, p, s))

  g0 <- proc.time()[["elapsed"]]
  x <- generate_data(n, p, s)$x
  gen_sec <- proc.time()[["elapsed"]] - g0

  error <- ""; res <- NULL; console <- character(0)
  t0 <- proc.time()[["elapsed"]]
  tryCatch({
    if (cfg$family == "TCV") {
      ## multiDT() reports progress with cat(); capture it, as run_icdt(quiet = TRUE)
      ## does for fac_ic(). Warnings are muffled for both methods alike.
      console <- utils::capture.output(
        res <- suppressWarnings(multiDT(x, K = cfg$K, rmax = RMAX)),
        type = "output")
    } else {
      ## Same call as in the ICDT comparison scripts (run_icdt_<scenario>.R).
      res <- suppressWarnings(run_icdt(x, qmax = RMAX, reptimes = cfg$M,
                                       maxit = MAXITER, quiet = TRUE,
                                       restore_rng = FALSE))
    }
  }, error = function(e) error <<- conditionMessage(e))
  elapsed <- proc.time()[["elapsed"]] - t0

  ## A call that returns is "ok" only if its criterion values are complete and
  ## finite (same rule as the ICDT comparison scripts).
  valid <- if (nzchar(error)) FALSE else if (cfg$family == "TCV") {
    is.list(res) && length(res$TCV) == RMAX && all(is.finite(res$TCV))
  } else {
    is.list(res) && length(res$values$ICDT_val) == RMAX &&
      all(is.finite(res$values$ICDT_val)) && length(res$ICDT) == 1L &&
      res$ICDT %in% seq_len(RMAX)
  }
  status <- if (nzchar(error)) "error" else if (valid) "ok" else "nonfinite_criterion"

  ## multiDT() itself prints "jic" every time a GFM fit fails and CJMLE is used
  ## instead; that count is read off its own output (NA if the call did not
  ## return, since the output is then lost). fac_ic() prints nothing when it
  ## switches to its VEM fallback, so for ICDT this is not available.
  fallback <- if (cfg$family == "TCV" && !nzchar(error)) sum(trimws(console) == "jic") else NA_integer_

  row <- list(method = cfg$method, family = cfg$family, K = cfg$K, M = cfg$M,
              n = n, p = p, s = s, seed = s + DATA_SEED_OFFSET, true_q = TRUE_Q,
              rmax = RMAX, maxiter = MAXITER,
              status = status, error = error,
              gen_sec = gen_sec,
              elapsed_sec = elapsed,                    # whole method call
              n_fits_nominal = tc_n_fits_nominal(cfg),  # K*8 or (M+1)*8, by construction
              fallback_fits = fallback,                 # TCV: GFM fits replaced by CJMLE; ICDT: NA
              started = started, finished = Sys.time(),
              host = Sys.info()[["nodename"]], pid = Sys.getpid(),
              workers = WORKERS, R_version = R.version.string,
              GFM_version = as.character(utils::packageVersion("GFM")))
  if (!nzchar(error) && is.list(res)) {
    if (cfg$family == "TCV") {
      row$curve_TCV <- res$TCV
      if (valid) {
        ## Same bookkeeping as the original simulation scripts: the four curves
        ## are built from re_TCV$TCV / re_TCV$TICV and which.min() picks the rank.
        re_TCV <- res
        log_TCV   <- re_TCV$TCV
        log_TPCV  <- re_TCV$TCV + seq(1:8) * ((n+p)/(n*p)) * log((n*p)/(n+p))
        log_TICV  <- re_TCV$TICV
        log_TIPCV <- re_TCV$TICV + seq(1:8) * ((n+p)/(n*p)) * log((n*p)/(n+p))
        row$sel_TCV   <- which.min(log_TCV)
        row$sel_TPCV  <- which.min(log_TPCV)
        row$sel_TICV  <- which.min(log_TICV)
        row$sel_TIPCV <- which.min(log_TIPCV)
      }
    } else {
      row$curve_ICDT <- res$values$ICDT_val
      row$fac_ic_time_sec <- res$time_sec     # fac_ic() call only, as timed inside run_icdt()
      if (valid) {
        row$sel_ICDT <- res$ICDT
        row$sel_VSDT <- res$VSDT
        row$sel_IC4  <- res$IC4
        row$sel_IC5  <- res$IC5
        row$VSDT_each <- res$VSDT_each
      }
    }
  }
  tc_save_row(row, path)
  tc_log(paths, sprintf("DONE  %-6s n=%4d p=%4d s=%3d %s seconds=%.1f fallback_fits=%s%s",
                        cfg$method, n, p, s, status, elapsed,
                        if (is.na(fallback)) "NA" else as.character(fallback),
                        if (nzchar(error)) paste0(" error=", error) else ""))
  tc_config_progress(paths, cfg, n, p)
  invisible(path)
}

## ---------------------------------------------------------------------------
## Collect rows
## ---------------------------------------------------------------------------
.tc_or_na <- function(v) if (is.null(v) || !length(v)) NA else v[[1]]

## Read result files into one data frame. Unreadable files are reported in the
## "unreadable" attribute and skipped, never silently dropped.
tc_read_rows <- function(files) {
  if (!length(files)) return(NULL)
  bad <- character(0)
  rows <- lapply(files, function(f) {
    z <- tryCatch(readRDS(f), error = function(e) NULL)
    if (is.null(z) || !is.list(z)) { bad <<- c(bad, f); return(NULL) }
    data.frame(method = z$method, family = z$family, K = .tc_or_na(z$K), M = .tc_or_na(z$M),
               n = z$n, p = z$p, s = z$s, seed = z$seed, status = z$status,
               gen_sec = z$gen_sec, elapsed_sec = z$elapsed_sec,
               n_fits_nominal = .tc_or_na(z$n_fits_nominal),
               fallback_fits = .tc_or_na(z$fallback_fits),
               sel_TCV = .tc_or_na(z$sel_TCV), sel_TPCV = .tc_or_na(z$sel_TPCV),
               sel_TICV = .tc_or_na(z$sel_TICV), sel_TIPCV = .tc_or_na(z$sel_TIPCV),
               sel_ICDT = .tc_or_na(z$sel_ICDT), sel_VSDT = .tc_or_na(z$sel_VSDT),
               sel_IC4 = .tc_or_na(z$sel_IC4), sel_IC5 = .tc_or_na(z$sel_IC5),
               host = z$host, workers = .tc_or_na(z$workers),
               started = format(z$started, "%Y-%m-%d %H:%M:%S"),
               error = z$error, folder = basename(dirname(dirname(f))),
               file = f, stringsAsFactors = FALSE)
  })
  rows <- rows[!vapply(rows, is.null, logical(1))]
  if (!length(rows)) { df <- NULL } else {
    df <- do.call(rbind, rows)
    df$method <- factor(df$method, levels = METHOD_ORDER)
    df$config <- sprintf("%dx%d", df$n, df$p)
    df <- df[order(df$method, -df$n * df$p, df$n, df$s), ]
  }
  attr(df, "unreadable") <- bad
  df
}

## rows.csv for one method, written by that method's own run only.
tc_collect_method <- function(paths) {
  df <- tc_read_rows(Sys.glob(file.path(paths$rows, "*.rds")))
  if (is.null(df)) return(invisible(NULL))
  write.csv(df, file.path(paths$method_dir, "rows.csv"), row.names = FALSE)
  invisible(df)
}

## ---------------------------------------------------------------------------
## Cross-method summary (summarize_time.R, after all six runs have finished)
## ---------------------------------------------------------------------------
tc_summarize <- function(root) {
  res_dir <- file.path(root, "results")
  df <- tc_read_rows(Sys.glob(file.path(res_dir, "*", "rows", "*.rds")))
  unreadable <- attr(df, "unreadable")
  if (is.null(df)) { message("No readable result files under ", res_dir); return(invisible(NULL)) }
  ## A file whose recorded method differs from the folder it sits in is reported
  ## and left out (a run would have stopped on it; summarize must not count it).
  misplaced <- df$file[df$folder != as.character(df$method)]
  df <- df[df$folder == as.character(df$method), ]
  if (!nrow(df)) { message("No result files in their own method folder under ", res_dir); return(invisible(NULL)) }
  methods <- METHOD_ORDER[METHOD_ORDER %in% as.character(unique(df$method))]
  configs <- sprintf("%dx%d", SIZES$n, SIZES$p)

  ## One row per method x designed configuration, built from ALL records, so
  ## that groups with no successful task are kept and missing tasks are counted.
  summ <- do.call(rbind, lapply(methods, function(m) do.call(rbind, lapply(seq_len(nrow(SIZES)), function(i) {
    d <- df[as.character(df$method) == m & df$n == SIZES$n[i] & df$p == SIZES$p[i], ]
    ok <- d[d$status == "ok", ]
    fb <- ok$fallback_fits
    data.frame(method = m, config = configs[i], n = SIZES$n[i], p = SIZES$p[i],
               designed = S, done = nrow(d), ok = nrow(ok),
               errors = sum(d$status == "error"),
               nonfinite = sum(d$status == "nonfinite_criterion"),
               missing = S - nrow(d),
               n_fits_nominal = if (nrow(d)) d$n_fits_nominal[1] else NA_integer_,
               median_sec = if (nrow(ok)) median(ok$elapsed_sec) else NA_real_,
               mean_sec = if (nrow(ok)) mean(ok$elapsed_sec) else NA_real_,
               sd_sec = if (nrow(ok) > 1) sd(ok$elapsed_sec) else NA_real_,
               min_sec = if (nrow(ok)) min(ok$elapsed_sec) else NA_real_,
               max_sec = if (nrow(ok)) max(ok$elapsed_sec) else NA_real_,
               tasks_with_fallback = if (!nrow(ok) || all(is.na(fb))) NA_integer_ else sum(fb > 0, na.rm = TRUE),
               fallback_fits_total = if (!nrow(ok) || all(is.na(fb))) NA_integer_ else sum(fb, na.rm = TRUE),
               median_gen_sec = if (nrow(d)) median(d$gen_sec) else NA_real_,
               workers = paste(sort(unique(d$workers)), collapse = ";"),
               hosts = paste(unique(d$host), collapse = ";"),
               stringsAsFactors = FALSE)
  }))))
  ## ratio of median time to TCV2 in the same configuration
  ref <- summ[summ$method == "TCV2", c("config", "median_sec")]
  names(ref)[2] <- "tcv2_median_sec"
  summ <- merge(summ, ref, by = "config", all.x = TRUE, sort = FALSE)
  summ$ratio_to_TCV2 <- summ$median_sec / summ$tcv2_median_sec
  summ <- summ[order(match(summ$config, configs), match(summ$method, METHOD_ORDER)), ]

  dir.create(file.path(res_dir, "summary"), recursive = TRUE, showWarnings = FALSE)
  write.csv(df, file.path(res_dir, "summary", "all_rows.csv"), row.names = FALSE)
  write.csv(summ, file.path(res_dir, "summary", "summary_by_config.csv"), row.names = FALSE)

  ## Markdown report: configurations x methods
  cell <- function(cfg, m, col, fmt) {
    v <- summ[summ$config == cfg & summ$method == m, col]
    if (!length(v) || is.na(v)) "" else sprintf(fmt, v)
  }
  tab <- function(col, fmt, title) {
    c(paste0("### ", title), "",
      paste0("| config | ", paste(methods, collapse = " | "), " |"),
      paste0("|---|", paste(rep("---:", length(methods)), collapse = "|"), "|"),
      vapply(configs, function(cfg) paste0("| ", cfg, " | ",
        paste(vapply(methods, function(m) cell(cfg, m, col, fmt), ""), collapse = " | "), " |"), ""),
      "")
  }
  completion <- vapply(configs, function(cfg) paste0("| ", cfg, " | ",
    paste(vapply(methods, function(m) {
      r <- summ[summ$config == cfg & summ$method == m, ]
      sprintf("%d/%d ok, %d not ok, %d missing", r$ok, r$designed, r$errors + r$nonfinite, r$missing)
    }, ""), collapse = " | "), " |"), "")
  fits <- vapply(methods, function(m) { v <- summ$n_fits_nominal[summ$method == m]; v[!is.na(v)][1] }, 0)
  hosts <- vapply(methods, function(m) paste(unique(unlist(strsplit(summ$hosts[summ$method == m], ";"))), collapse = ";"), "")
  workers <- vapply(methods, function(m) paste(unique(unlist(strsplit(summ$workers[summ$method == m], ";"))), collapse = ";"), "")
  report <- c("# TCV vs ICDT computing time (Example 1)", "",
    sprintf("Generated %s from %s.", format(Sys.time(), "%Y-%m-%d %H:%M"), res_dir), "",
    "Time = the whole method call for one count matrix (seconds, single-threaded process): thinning, all model fits",
    "including any fallback the method performs itself, predicted intensities and held-out losses. Data generation excluded.",
    "Both methods use GFM::gfm (Poisson, AM, maxIter 30) as their main solver and keep their own fallback",
    "(TCV: CJMLE after a failed GFM fit; ICDT: the authors' VEM re-run of a block of 8 fits after a GFM error).",
    "'Fits' is the nominal number implied by the method's structure (K x 8 for TCV_K, (M+1) x 8 for ICDT_M), not a count of solver calls.", "",
    "| method | nominal fits per data set | workers | hosts |", "|---|---:|---|---|",
    sprintf("| %s | %s | %s | %s |", methods, ifelse(is.na(fits), "", fits), workers, hosts), "",
    "### Completion", "",
    paste0("| config | ", paste(methods, collapse = " | "), " |"),
    paste0("|---|", paste(rep("---", length(methods)), collapse = "|"), "|"),
    completion, "",
    tab("median_sec", "%.1f", "Median seconds per data set (ok tasks)"),
    tab("ratio_to_TCV2", "%.2f", "Ratio of median seconds to TCV2"),
    tab("tasks_with_fallback", "%d", "Tasks with at least one fallback fit (TCV only; not observable for ICDT)"),
    if (length(unreadable)) c("### Unreadable result files (skipped)", "", paste0("- ", unreadable), "") else character(0),
    if (length(misplaced)) c("### Result files in the wrong method folder (skipped)", "", paste0("- ", misplaced), "") else character(0),
    "Rows with status error or nonfinite_criterion are excluded from the time statistics; missing = designed tasks without a result file.")
  writeLines(report, file.path(res_dir, "summary", "REPORT.md"))
  cat("\nSUMMARY written to", file.path(res_dir, "summary"), "\n")
  if (length(unreadable)) cat("Unreadable result files skipped:", length(unreadable), "\n")
  if (length(misplaced)) cat("Result files in the wrong method folder skipped:", length(misplaced), "\n")
  print(summ[, c("method", "config", "ok", "errors", "nonfinite", "missing", "n_fits_nominal",
                 "median_sec", "ratio_to_TCV2", "tasks_with_fallback")], row.names = FALSE, digits = 4)
  invisible(summ)
}

## ---------------------------------------------------------------------------
## Entry point used by the six run scripts
## ---------------------------------------------------------------------------
tc_main <- function(cfg) {
  paths <- tc_paths(cfg$root, cfg$method)
  dir.create(paths$rows, recursive = TRUE, showWarnings = FALSE)
  f <- tc_files(cfg$root)
  if (cfg$family == "TCV") {
    tc_check_packages(f$core)
    tc_load_tcv(f)
  } else {
    tc_check_packages(c(f$core, f$wrapper, f$functions))
    tc_load_icdt(f, cfg$root)
  }
  tc_load_generator(f)
  capture.output(sessionInfo(), file = file.path(paths$method_dir, "sessionInfo.txt"))

  tasks <- merge(SIZES, data.frame(s = seq_len(S)))
  ## Same order for every method: the order of SIZES (unbalanced configurations
  ## last), then replicate.
  tasks <- tasks[order(match(paste(tasks$n, tasks$p), paste(SIZES$n, SIZES$p)), tasks$s), ]
  task_paths <- vapply(seq_len(nrow(tasks)), function(i)
    tc_task_path(paths, tasks$n[i], tasks$p[i], tasks$s[i]), "")

  ## Leftovers of interrupted writes are never read as results.
  tmp <- list.files(paths$rows, pattern = "\\.rds\\.tmp$", full.names = TRUE)
  if (length(tmp))
    tc_log(paths, sprintf("NOTE %d interrupted temporary file(s) ignored: %s",
                          length(tmp), paste(basename(tmp), collapse = " ")))
  ## Existing results are accepted only after validation (tc_check_existing stops otherwise).
  exists <- file.exists(task_paths)
  same_workers <- vapply(which(exists), function(i)
    tc_check_existing(task_paths[i], cfg, tasks$n[i], tasks$p[i], tasks$s[i]), logical(1))
  if (any(!same_workers))
    tc_log(paths, sprintf("WARNING %d existing result(s) were produced with a different TIME_WORKERS than the current %d; say so when comparing times.",
                          sum(!same_workers), WORKERS))
  tc_log(paths, sprintf("%s tasks=%d already_done=%d workers=%d host=%s %s",
                        cfg$method, nrow(tasks), sum(exists), WORKERS,
                        Sys.info()[["nodename"]], R.version.string))

  todo <- tasks[!exists, , drop = FALSE]
  if (nrow(todo)) {
    task_list <- split(todo, seq_len(nrow(todo)))
    if (WORKERS > 1L) {
      out <- parallel::mclapply(task_list, function(t)
        tryCatch(tc_run_task(t, cfg, paths), error = function(e) e),
        mc.cores = WORKERS, mc.preschedule = FALSE)
      msg_of <- function(o) {
        if (is.null(o)) "child process returned NULL (crashed?)"
        else if (inherits(o, "try-error")) trimws(as.character(o))
        else if (inherits(o, "condition")) conditionMessage(o)
        else NA_character_
      }
      msgs <- vapply(out, msg_of, "")
      bad <- !is.na(msgs)
      if (any(bad))
        tc_log(paths, sprintf("WARNING %d task(s) left no result file. First message: %s",
                              sum(bad), msgs[bad][1]))
    } else {
      for (t in task_list) tc_run_task(t, cfg, paths)
    }
  }

  ## Compare with the design: every (n, p, s) must have a result file now.
  missing <- !file.exists(task_paths)
  tc_collect_method(paths)
  tc_log(paths, sprintf("%s finished: %d of %d tasks have result files; rows in %s; rows.csv in %s. Run summarize_time.R after all six scripts have finished.",
                        cfg$method, sum(!missing), length(task_paths), paths$rows, paths$method_dir))
  if (any(missing)) {
    shown <- head(basename(task_paths[missing]), 10)
    tc_log(paths, sprintf("MISSING %d task(s) without result file, e.g. %s. Rerun the same command to retry them.",
                          sum(missing), paste(shown, collapse = " ")))
    stop(sum(missing), " task(s) of ", cfg$method, " have no result file", call. = FALSE)
  }
  invisible(TRUE)
}
