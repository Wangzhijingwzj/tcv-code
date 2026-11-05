## make_tables.R -- summary tables for the timing study.
## Computes every statistic from the raw per-task files results/<METHOD>/rows/*.rds,
## cross-checks them against results/summary/{all_rows.csv, summary_by_config.csv, REPORT.md},
## and writes output/*.csv, output/*.tex and output/checks.txt. Read-only on results/.
## Any failed check stops the script before tables are written.

tc_script_dir <- function() {
  f <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE))
  if (length(f)) dirname(normalizePath(f[1])) else normalizePath(getwd())
}
ROOT <- tc_script_dir()
RES  <- file.path(ROOT, "results")
OUT  <- file.path(ROOT, "output")
dir.create(OUT, showWarnings = FALSE)

METHODS <- c("TCV2", "TCV5", "TCV10", "ICDT3", "ICDT5", "ICDT10")
LABELS  <- c("TCV K=2", "TCV K=5", "TCV K=10", "ICDT M=3", "ICDT M=5", "ICDT M=10")
NFITS   <- c(TCV2 = 16L, TCV5 = 40L, TCV10 = 80L, ICDT3 = 32L, ICDT5 = 48L, ICDT10 = 88L)
SIZES   <- data.frame(n = c(100L, 200L, 500L, 1000L, 100L, 1000L),
                      p = c(100L, 200L, 500L, 1000L, 1000L, 100L))
CONFIGS <- sprintf("%dx%d", SIZES$n, SIZES$p)
S <- 100L

checks <- character(0)
check <- function(ok, msg) {
  line <- sprintf("%s  %s", if (isTRUE(ok)) "PASS" else "FAIL", msg)
  checks <<- c(checks, line)
  cat(line, "\n")
  if (!isTRUE(ok)) { writeLines(checks, file.path(OUT, "checks.txt")); stop("check failed: ", msg) }
}

## 1. raw per-task result files -------------------------------------------------
files <- unlist(lapply(METHODS, function(m) Sys.glob(file.path(RES, m, "rows", "*.rds"))))
n_expected <- length(METHODS) * nrow(SIZES) * S
check(length(files) == n_expected, sprintf("%d result files found (expected %d)", length(files), n_expected))
fields <- c("method", "n", "p", "s", "seed", "status", "elapsed_sec", "gen_sec", "n_fits_nominal", "fallback_fits")
raw <- do.call(rbind, lapply(files, function(f) {
  x <- readRDS(f)
  v <- lapply(fields, function(k) if (is.null(x[[k]])) NA else x[[k]])
  names(v) <- fields
  as.data.frame(v, stringsAsFactors = FALSE)
}))
raw$config <- sprintf("%dx%d", raw$n, raw$p)
check(all(raw$status == "ok"), "every task has status ok")
check(all(raw$method %in% METHODS) && all(raw$config %in% CONFIGS), "only the designed methods and configurations occur")
cnt <- table(factor(raw$method, METHODS), factor(raw$config, CONFIGS))
check(all(cnt == S), sprintf("every method x configuration has exactly %d tasks", S))
same_s <- tapply(raw$s, list(raw$method, raw$config), function(v) identical(sort(as.integer(v)), seq_len(S)))
check(all(same_s), "within every cell the replicate index s is exactly 1..100 (same data sets for all methods)")
check(all(raw$seed == 1000L + raw$s), "seed == 1000 + s for every task")
check(all(raw$n_fits_nominal == NFITS[raw$method]), "nominal fit counts are 8K (TCV) and 8(M+1) (ICDT)")
check(all(is.finite(raw$elapsed_sec)) && all(raw$elapsed_sec > 0), "elapsed_sec is finite and positive everywhere")

## 2. cross-check against all_rows.csv ---------------------------------------------
ar <- read.csv(file.path(RES, "summary", "all_rows.csv"), stringsAsFactors = FALSE)
key <- function(d) paste(d$method, d$n, d$p, d$s)
m <- match(key(raw), key(ar))
check(nrow(ar) == nrow(raw) && !anyNA(m) && !any(duplicated(m)), "all_rows.csv contains exactly the same 3600 tasks")
check(max(abs(raw$elapsed_sec - ar$elapsed_sec[m])) < 1e-9 && all(raw$status == ar$status[m]),
      "elapsed_sec and status in all_rows.csv equal the raw files")

## 3. statistics per method x configuration ----------------------------------------
cells <- expand.grid(config = CONFIGS, method = METHODS, stringsAsFactors = FALSE)
stat <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  d <- raw[raw$method == cells$method[i] & raw$config == cells$config[i], ]
  x <- d$elapsed_sec
  fb <- d$fallback_fits
  data.frame(method = cells$method[i], config = cells$config[i], n = d$n[1], p = d$p[1],
             n_fits_nominal = NFITS[[cells$method[i]]], n_ok = length(x),
             mean_sec = mean(x), sd_sec = sd(x),
             median_sec = median(x),
             mad_sec = mad(x, constant = 1),      # median(|x - median(x)|), unscaled
             mad_scaled_sec = mad(x),             # x 1.4826, comparable with sd under normality
             q1_sec = unname(quantile(x, 0.25)), q3_sec = unname(quantile(x, 0.75)),
             min_sec = min(x), max_sec = max(x),
             median_sec_per_fit = median(x) / NFITS[[cells$method[i]]],
             tasks_with_fallback = if (all(is.na(fb))) NA_integer_ else sum(fb > 0, na.rm = TRUE),
             fallback_fits_total = if (all(is.na(fb))) NA_integer_ else sum(fb, na.rm = TRUE),
             stringsAsFactors = FALSE)
}))
ref <- stat[stat$method == "TCV2", ]
stat$ratio_median_to_TCV2 <- stat$median_sec / ref$median_sec[match(stat$config, ref$config)]
stat$ratio_mean_to_TCV2   <- stat$mean_sec   / ref$mean_sec[match(stat$config, ref$config)]
stat <- stat[order(match(stat$config, CONFIGS), match(stat$method, METHODS)), ]
rownames(stat) <- NULL

## 4. cross-check against summary_by_config.csv and REPORT.md ------------------------
sb <- read.csv(file.path(RES, "summary", "summary_by_config.csv"), stringsAsFactors = FALSE)
j <- match(paste(stat$method, stat$config), paste(sb$method, sb$config))
check(nrow(sb) == nrow(stat) && !anyNA(j), "summary_by_config.csv has the same 36 cells")
check(max(abs(stat$median_sec - sb$median_sec[j])) < 1e-6 &&
      max(abs(stat$mean_sec - sb$mean_sec[j])) < 1e-6 &&
      max(abs(stat$sd_sec - sb$sd_sec[j])) < 1e-6 &&
      all(stat$n_ok == sb$ok[j]) &&
      identical(as.integer(stat$tasks_with_fallback), as.integer(sb$tasks_with_fallback[j])),
      "median, mean, sd, ok counts and fallback counts equal summary_by_config.csv")

rep <- readLines(file.path(RES, "summary", "REPORT.md"))
start <- grep("^### Median seconds", rep)
hdr <- trimws(strsplit(rep[start + 2], "\\|")[[1]])[-1]
rows <- rep[(start + 4):(start + 3 + length(CONFIGS))]
ok_rep <- length(start) == 1 && all(sapply(rows, function(l) {
  r <- trimws(strsplit(l, "\\|")[[1]])[-1]
  all(sapply(2:length(hdr), function(k)
    r[k] == sprintf("%.1f", stat$median_sec[stat$config == r[1] & stat$method == hdr[k]])))
}))
check(ok_rep, "medians agree with the median table printed in REPORT.md (1 decimal)")

## 5. write outputs ------------------------------------------------------------------
num_cols <- sapply(stat, is.numeric)
stat_out <- stat; stat_out[num_cols] <- lapply(stat_out[num_cols], round, 4)
write.csv(stat_out, file.path(OUT, "stats_by_config.csv"), row.names = FALSE)

cell_tab <- function(a, b, fmt = "%.2f (%.2f)") {
  w <- sapply(METHODS, function(mm) sapply(CONFIGS, function(cf) {
    r <- stat[stat$method == mm & stat$config == cf, ]; sprintf(fmt, r[[a]], r[[b]]) }))
  w <- as.data.frame(w, stringsAsFactors = FALSE); names(w) <- LABELS
  cbind(config = sprintf("(%d, %d)", SIZES$n, SIZES$p), w, stringsAsFactors = FALSE)
}
num_tab <- function(a) {
  w <- sapply(METHODS, function(mm) sapply(CONFIGS, function(cf) stat[[a]][stat$method == mm & stat$config == cf]))
  w <- as.data.frame(round(w, 2)); names(w) <- LABELS
  cbind(config = sprintf("(%d, %d)", SIZES$n, SIZES$p), w, stringsAsFactors = FALSE)
}
tab_mm <- cell_tab("median_sec", "mad_sec")
tab_ms <- cell_tab("mean_sec", "sd_sec")
write.csv(tab_mm, file.path(OUT, "table_median_mad.csv"), row.names = FALSE)
write.csv(tab_ms, file.path(OUT, "table_mean_sd.csv"),   row.names = FALSE)
write.csv(num_tab("median_sec"), file.path(OUT, "median_wide.csv"), row.names = FALSE)
write.csv(num_tab("mad_sec"),    file.path(OUT, "mad_wide.csv"),    row.names = FALSE)
write.csv(num_tab("mean_sec"),   file.path(OUT, "mean_wide.csv"),   row.names = FALSE)
write.csv(num_tab("sd_sec"),     file.path(OUT, "sd_wide.csv"),     row.names = FALSE)

tex_table <- function(tab, stat_label, label, file) {
  body <- apply(tab, 1, function(r) sprintf("$%s$ & %s \\\\", gsub(" ", "", r[1]), paste(r[-1], collapse = " & ")))
  lines <- c(
    "\\begin{table}[htbp]", "\\centering",
    sprintf(paste0("\\caption{Computing time in seconds for one $n\\times p$ count matrix in Example~1: ",
                   "%s over 100 replicates. TCV with $K$ folds fits $8K$ GFM models and ICDT with $M$ thinning ",
                   "repetitions fits $8(M+1)$ models ($r_{\\max}=q_{\\max}=8$); both methods use \\texttt{GFM::gfm} ",
                   "(Poisson, AM algorithm, 30 iterations) as the main solver, and each method's own fallback is ",
                   "included in the time. Every replicate ran as a single-threaded process.}"), stat_label),
    sprintf("\\label{%s}", label),
    "\\footnotesize", "\\setlength{\\tabcolsep}{3pt}",
    "\\begin{tabular}{lrrrrrr}", "\\toprule",
    " & \\multicolumn{3}{c}{TCV ($K$ folds)} & \\multicolumn{3}{c}{ICDT ($M$ thinnings)} \\\\",
    "\\cmidrule(lr){2-4} \\cmidrule(lr){5-7}",
    "$(n,p)$ & $K=2$ & $K=5$ & $K=10$ & $M=3$ & $M=5$ & $M=10$ \\\\",
    "\\midrule", body[1:4], "\\addlinespace", body[5:6],
    "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(lines, file)
}
tex_table(tab_mm, "median (median absolute deviation)", "tab:time-median", file.path(OUT, "table_median_mad.tex"))
tex_table(tab_ms, "mean (standard deviation)",          "tab:time-mean",   file.path(OUT, "table_mean_sd.tex"))

writeLines(c(sprintf("make_tables.R run %s on %d result files", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), nrow(raw)), checks,
             "", "MAD = median(|x - median(x)|), unscaled (R: mad(x, constant = 1)); mad_scaled_sec in stats_by_config.csv uses the default factor 1.4826."),
           file.path(OUT, "checks.txt"))

cat("\nAll checks passed. Files written to", OUT, ":\n")
cat(paste0("  ", list.files(OUT)), sep = "\n")
cat("\nTable: median (MAD), seconds\n")
print(tab_mm, row.names = FALSE)
