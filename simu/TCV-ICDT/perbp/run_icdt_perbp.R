# ICDT comparison: original TCV data generation; data seeds 1000 + seq_len(S).
# Run: Rscript --vanilla run_icdt_perbp.R
# Resume: run the same command; completed replicate files are retained.
# Other modes: summarize (saved results only).
ROOT <- path.expand("~/multiDT/TCV-ICDT")
SCENARIO <- "perbp"
SOURCE_REL <- "perbp/source/250410_TCV_perbp.R"
HISTORY_PATTERN <- "fnorm_perbp_n%d_p%d_q%d_with_counts.rds"
TRUE_Q <- 5L
WORK_DIR <- file.path(ROOT, SCENARIO)
SCRIPT_PATH <- file.path(WORK_DIR, paste0("run_icdt_", SCENARIO, ".R"))
RESULT_DIR <- file.path(WORK_DIR, "results")
REPLICATE_DIR <- file.path(RESULT_DIR, "replicates")
WRAPPER <- file.path(ROOT, "common/ICDT/icdt_wrapper.R")
CORE <- file.path(ROOT, "common/250113function_tcv.R")
S <- 100L
M <- 3L
QMAX <- 8L
MAXITER <- 30L
WORKERS <- 39L
SIZES <- c(100L, 200L, 500L, 1000L)
Sys.setenv(OMP_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1",
           MKL_NUM_THREADS = "1", VECLIB_MAXIMUM_THREADS = "1")
suppressPackageStartupMessages({ library(irlba); library(Rfast) })

# Read only the original generator block, without running its simulation loop.
# Three separate add_identifiability calls, their RNG consumption, the Poisson
# loops and (where present) the zero-inflation rejection loop are kept intact.
source_lines <- readLines(file.path(ROOT, SOURCE_REL), warn = FALSE)
first <- grep("^[[:space:]]*set.seed\\(s \\+ 1000\\)", source_lines)
last <- grep("^[[:space:]]*all_results\\$x <- x", source_lines)
GENERATION_EXPR <- parse(text = c("all_results <- list()", source_lines[first:last]))
helper <- Filter(function(e) is.call(e) &&
                   as.character(e[[1]]) %in% c("<-", "=") &&
                   identical(e[[2]], as.name("add_identifiability")),
                 as.list(parse(CORE)))
GEN_ENV <- new.env(parent = globalenv())
eval(helper[[1]], envir = GEN_ENV)
GEN_ENV$Diag <- pracma::Diag
GEN_ENV$irlba <- irlba::irlba
GEN_ENV$has_full_zero_row_or_col <- function(mat)
  any(rowSums(mat == 0) == ncol(mat)) || any(colSums(mat == 0) == nrow(mat))

generate_data <- function(n, p, s) {
  e <- new.env(parent = GEN_ENV)
  e$n <- n; e$p <- p; e$s <- s
  eval(GENERATION_EXPR, envir = e)
  list(x = e$x, q = e$q)
}

checkpoint_path <- function(n, p, s)
  file.path(REPLICATE_DIR, sprintf("n%d_p%d_seed%d.rds", n, p, 1000L + s))

init_icdt <- function() {
  if (!exists("run_icdt", envir = globalenv(), inherits = FALSE)) {
    source(WRAPPER, local = globalenv())
  }
}

run_one <- function(task) {
  n <- as.integer(task$n); p <- as.integer(task$p); s <- as.integer(task$s)
  path <- checkpoint_path(n, p, s)
  if (file.exists(path)) {
    old <- readRDS(path)
    stopifnot(old$n == n, old$p == p, old$seed == 1000L+s,
              old$M == M, old$qmax == QMAX, old$maxIter == MAXITER,
              old$scenario == SCENARIO, old$q == TRUE_Q)
    return(path)
  }
  started <- Sys.time()
  cat(format(started), SCENARIO, "START", n, p, 1000L+s, "\n")
  x <- NULL; result <- NULL; error <- ""
  tryCatch(suppressWarnings({
    x <- generate_data(n, p, s)$x
    result <- run_icdt(x, qmax = QMAX, reptimes = M, maxit = MAXITER,
                       quiet = TRUE, restore_rng = FALSE)
  }), error = function(e) error <<- conditionMessage(e))
  valid <- !is.null(result) && length(result$values$ICDT_val) == QMAX &&
    all(is.finite(result$values$ICDT_val)) && length(result$ICDT) == 1L &&
    result$ICDT %in% seq_len(QMAX)
  status <- if (nzchar(error)) "error" else if (valid) "ok" else "nonfinite_criterion"
  record <- list(scenario = SCENARIO, n = n, p = p, s = s, seed = 1000L+s,
    q = TRUE_Q, M = M, qmax = QMAX, maxIter = MAXITER, status = status,
    started = started, finished = Sys.time(),
    elapsed_sec = as.numeric(difftime(Sys.time(), started, units = "secs")),
    error = error, x = x, result = result)
  saveRDS(record, path)
  cat(format(Sys.time()), SCENARIO, "DONE", n, p, 1000L+s, status,
      "ICDT", if (valid) result$ICDT else NA_integer_,
      "seconds", round(record$elapsed_sec, 1),
      "\n")
  path
}

stat_row <- function(n, p, method, values, attempted, source) {
  valid <- is.finite(values) & values >= 1 & values <= QMAX
  v <- floor(values[valid])
  co <- sum(v == TRUE_Q); ov <- sum(v > TRUE_Q); un <- sum(v < TRUE_Q)
  data.frame(n = n, p = p, method = method, true_q = TRUE_Q, requested = S,
    valid = length(v), failed = attempted-length(v),
    correct = co, over = ov, under = un,
    MO = if (ov) mean(v[v > TRUE_Q]-TRUE_Q) else 0,
    MU = if (un) mean(TRUE_Q-v[v < TRUE_Q]) else 0,
    source = source, stringsAsFactors = FALSE)
}

summarize_results <- function() {
  grid <- expand.grid(n = SIZES, p = SIZES, s = seq_len(S))
  records <- list(); selected <- list()
  for (i in seq_len(nrow(grid))) {
    task <- grid[i, ]; path <- checkpoint_path(task$n, task$p, task$s)
    if (!file.exists(path)) next
    z <- readRDS(path)
    r <- z$result
    pick <- function(name, curve) {
      if (is.null(r) || length(r[[name]]) != 1L ||
          any(!is.finite(r$values[[curve]]))) return(NA_integer_)
      as.integer(r[[name]])
    }
    row <- data.frame(n=z$n, p=z$p, s=z$s, seed=z$seed, status=z$status,
      ICDT=if (z$status == "ok") pick("ICDT", "ICDT_val") else NA_integer_,
      VSDT=pick("VSDT", "VSDT_val_mat"), IC4=pick("IC4", "IC4_val"),
      IC5=pick("IC5", "IC5_val"),
      error=z$error, stringsAsFactors=FALSE)
    records[[length(records)+1L]] <- row
  }
  if (!length(records)) { cat("No completed replicates yet.\n"); return(invisible(NULL)) }
  runs <- do.call(rbind, records)
  historical <- list(); summary <- list(); source_files <- character()
  mapping <- c(BCV="facnum_bcv", CV="facnum_bcv2", SR="facnum_sr",
    Ratio="facnum_ratio", JIC1="facnum_IC1", JIC2="facnum_IC2",
    TCV2="facnum_TCV2K", TPCV2="facnum_TPCV2K", TIPCV2="facnum_TIPCV2K")
  for (n in SIZES) for (p in SIZES) {
    rows <- runs[runs$n == n & runs$p == p, ]
    for (method in c("ICDT", "VSDT", "IC4", "IC5"))
      summary[[length(summary)+1L]] <- stat_row(n,p,method,rows[[method]],nrow(rows),"new_M3")
    fname <- sprintf(HISTORY_PATTERN,n,p,TRUE_Q)
    old_path <- file.path(WORK_DIR,"history",fname)
    if (file.exists(old_path)) {
      old <- readRDS(old_path)
      stopifnot(nrow(old) >= 8L+S)
      source_files <- c(source_files,old_path)
      for (method in names(mapping)) {
        values <- as.numeric(old[8L+seq_len(S),mapping[[method]]])
        historical[[length(historical)+1L]] <- data.frame(n=n,p=p,s=seq_len(S),
          seed=1000L+seq_len(S),method=method,selected=values,source=old_path)
        summary[[length(summary)+1L]] <- stat_row(n,p,method,values,S,paste0("historical_first",S))
      }
    } else {
      stop("Historical with_counts file not found: ", old_path)
    }
    add_path <- file.path(WORK_DIR,"history",paste0("add_methods_",fname))
    if (file.exists(add_path)) {
      old <- readRDS(add_path)
      if (nrow(old) >= 8L+S && "facnum_ne" %in% colnames(old)) {
        values <- as.numeric(old[8L+seq_len(S),"facnum_ne"])
        source_files <- c(source_files,add_path)
        historical[[length(historical)+1L]] <- data.frame(n=n,p=p,s=seq_len(S),
          seed=1000L+seq_len(S),method="SE",selected=values,source=add_path)
        summary[[length(summary)+1L]] <- stat_row(n,p,"SE",values,S,paste0("historical_first",S))
      } else if (nrow(old) < 8L+S) {
        stop("Historical add_methods_ file has fewer than ", 8L+S, " rows: ", add_path)
      } else {
        stop("Historical add_methods_ file lacks facnum_ne column: ", add_path)
      }
    } else {
      stop("Historical add_methods_ file not found: ", add_path)
    }
  }
  summary <- do.call(rbind,summary)
  icdt <- summary[summary$method == "ICDT", ]
  stamp <- format(Sys.time(), "%Y%m%d-%H%M%S")
  dest <- file.path(RESULT_DIR,paste0("summary_",stamp))
  stopifnot(!dir.exists(dest)); dir.create(dest,recursive=TRUE)
  write.csv(icdt,file.path(dest,"ICDT_summary.csv"),row.names=FALSE)
  write.csv(summary,file.path(dest,"comparison.csv"),row.names=FALSE)
  write.csv(runs,file.path(dest,"replicates.csv"),row.names=FALSE)
  notes <- c(paste("Scenario:",SCENARIO),paste("Original generator:",SOURCE_REL),
    sprintf("n,p: %s; %d replicates; data seeds 1001:%d; rmax=%d; M=%d; maxIter=%d.", paste(SIZES,collapse=","),S,1000L+S,QMAX,M,MAXITER),
    "ICDT uses the authors' fac_ic() unchanged, including its AM-to-VEM fallback.",
    sprintf("Thinning seeds are set internally by fac_ic() to 1001:%d, independently of TCV draws.",1000L+M),
    sprintf("Historical methods are read from rows 9:%d of the matching original with_counts files; they are not rerun.",8L+S),
    "MO/MU: mean over/under-selection magnitude conditional on that event; 0 if no such event.",
    "Finite results with numerical warnings are retained; success is not a proof of optimizer convergence.")
  openxlsx::write.xlsx(list(ICDT_summary=icdt,comparison=summary,
    new_replicates=runs,historical_replicates=if(length(historical)) do.call(rbind,historical) else data.frame(),
    notes=data.frame(note=notes),sources=data.frame(path=source_files)),
    file.path(dest,"RESULTS.xlsx"),overwrite=FALSE)
  capture.output(sessionInfo(),file=file.path(dest,"sessionInfo.txt"))
  report <- c(paste0("# ",SCENARIO,"：ICDT 比较"),"",
    sprintf("已完成 %d/%d 次；ICDT 有效 %d 次，失败或非有限 %d 次。",nrow(runs),nrow(grid),sum(is.finite(runs$ICDT)),sum(!is.finite(runs$ICDT))),"",
    sprintf("每组 %d 次，数据种子 1001–%d，ICDT 的 M=%d，候选因子数 1–%d。",S,1000L+S,M,QMAX), "",
    "| n | p | 正确 | 高估 | 低估 | 失败 | MO | MU |",
    "|---:|---:|---:|---:|---:|---:|---:|---:|")
  for(i in seq_len(nrow(icdt))) report <- c(report,sprintf("| %d | %d | %d | %d | %d | %d | %.2f | %.2f |",
    icdt$n[i],icdt$p[i],icdt$correct[i],icdt$over[i],icdt$under[i],icdt$failed[i],icdt$MO[i],icdt$MU[i]))
  report <- c(report,"",sprintf("`RESULTS.xlsx` 包含与历史方法相同种子编号前 %d 次的比较，以及逐次选秩和来源。",S),
    "历史方法仅从旧结果读取，本轮不重新拟合。MO/MU 沿用原表的事件条件平均误差定义，无事件时记 0。")
  writeLines(report,file.path(dest,"REPORT_CN.md"))
  cat("SUMMARY",dest,"\n"); print(icdt[,c("n","p","correct","over","under","failed","MO","MU")],row.names=FALSE)
  invisible(dest)
}

main <- function(mode="run") {
  dir.create(REPLICATE_DIR,recursive=TRUE,showWarnings=FALSE)
  if (mode == "summarize") return(summarize_results())
  stopifnot(mode == "run")
  # schedule the larger jobs first to avoid a
  # long single-worker tail. Scheduling does not alter per-replicate seeds.
  tasks <- expand.grid(n=SIZES,p=SIZES,s=seq_len(S))
  tasks <- tasks[order(-tasks$n*tasks$p,tasks$n,tasks$p,tasks$s),]
  done <- vapply(seq_len(nrow(tasks)),function(i)
    file.exists(checkpoint_path(tasks$n[i],tasks$p[i],tasks$s[i])),logical(1))
  tasks <- tasks[!done,]
  if (nrow(tasks)) {
    log <- file.path(RESULT_DIR,paste0("workers_",format(Sys.time(),"%Y%m%d-%H%M%S"),".log"))
    cl <- parallel::makePSOCKcluster(WORKERS,outfile=log)
    on.exit(parallel::stopCluster(cl),add=TRUE)
    parallel::clusterCall(cl,function(path) {
      source(path,local=globalenv()); init_icdt(); TRUE
    },SCRIPT_PATH)
    parallel::parLapplyLB(cl,split(tasks,seq_len(nrow(tasks))),run_one,chunk.size=1L)
  }
  summarize_results()
}

if (sys.nframe() == 0L) {
  args <- commandArgs(trailingOnly=TRUE)
  main(if(length(args)) args[[1]] else "run")
}
