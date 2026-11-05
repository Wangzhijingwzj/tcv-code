# Shared settings for the NB example. No original example is modified.
# Run every script from this folder (Slurm: submit from here; R: setwd() to this folder first).
work_dir <- getwd()
N_VALUES <- c(100L, 200L, 500L, 1000L)
P_VALUES <- c(100L, 200L, 500L, 1000L)
# Optional argument splitting the grid into two jobs, e.g. Rscript --vanilla 250410_TCV_nb.R A
#   A: p = 100, 200, 500;  B: p = 1000;  no argument: all p (used by out_all_2510.R).
part <- commandArgs(trailingOnly = TRUE)
if (identical(part, "A")) P_VALUES <- c(100L, 200L, 500L)
if (identical(part, "B")) P_VALUES <- 1000L
S <- 100L         # replicates per (n, p); data seeds 1001, ..., 1000 + S
WORKERS <- 39L    # parallel worker processes = cores requested in the Slurm file minus 1
result_dir <- file.path(work_dir, paste0("results_S", S))
KS <- 2L
TRUE_Q <- 5L
RMAX <- 8L
NB_SIZE <- 10
MAXITER <- 30L
ICDT_M <- 3L
stopifnot(S >= 1L, WORKERS >= 1L, all(N_VALUES > RMAX), all(P_VALUES > RMAX),
          identical(KS,2L))
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1")
dir.create(result_dir, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(result_dir, "replicates"), showWarnings=FALSE)
dir.create(file.path(result_dir, "raw_results"), showWarnings=FALSE)

run_spec <- list(scenario="nb", size=NB_SIZE, true_q=TRUE_Q, rmax=RMAX,
                 K=KS, maxIter=MAXITER, loss="Poisson_NLL_plus_1e-20",
                 pilot="full-data_JIC2_exact_NLL", size_bounds=c(1e-4, 1e6))
rep_path <- function(n, p, s, prefix="main")
  file.path(result_dir, "replicates", sprintf("%s_n%d_p%d_seed%d.rds", prefix,n,p,1000L+s))
case_stem <- function(n,p,prefix="") sprintf("%sfnorm_nb_n%d_p%d_q%d",prefix,n,p,TRUE_Q)
save_new <- function(object,path) {
  if (file.exists(path)) stop("Existing file retained: ",path)
  saveRDS(object,path)
}
rescore_validation_record <- function(z) {
  if(identical(z$spec$loss,"exact_Poisson_NLL")) {
    z$original_run_spec <- z$spec
    if(!is.null(z$thinning)) {
      z$selected_exact_previous <- z$selected
      for(nm in names(z$thinning)) {
        rr <- z$thinning[[nm]]
        stopifnot(is.matrix(rr$legacy_fold_loss),nrow(rr$legacy_fold_loss)==RMAX)
        cv <- rowSums(rr$legacy_fold_loss)
        stopifnot(isTRUE(all.equal(cv,rr$legacy_TCV,check.attributes=FALSE)))
        rr$exact_TCV <- rr$TCV
        rr$exact_fold_loss <- rr$fold_loss
        rr$TCV <- cv; rr$TICV <- log(cv); rr$fold_loss <- rr$legacy_fold_loss
        rr$loss <- "Poisson_NLL_plus_1e-20"
        z$thinning[[nm]] <- rr
        scheme <- sub("_K[0-9]+$","",nm)
        prefix <- if(scheme=="Poisson") "" else paste0(scheme,"_")
        pen <- seq_len(RMAX)*((z$n+z$p)/(z$n*z$p))*log(z$n*z$p/(z$n+z$p))
        curves <- list(TCV=cv,TPCV=cv+pen,TICV=log(cv),TIPCV=log(cv)+pen)
        for(metric in names(curves)) {
          z[[paste0("log_",prefix,metric,rr$K,"K")]] <- curves[[metric]]
          z$selected[paste0("facnum_",prefix,metric,rr$K,"K")] <- pick_rank(curves[[metric]])
        }
      }
      z$status <- if(anyNA(z$selected)) "partial_failure" else "completed"
    }
    z$spec$loss <- "Poisson_NLL_plus_1e-20"
    z$score_update <- "TCV scores recovered from saved +1e-20 losses; no data regeneration or refitting."
  }
  z
}
check_record <- function(z,n,p,s) {
  z <- rescore_validation_record(z)
  stopifnot(identical(z$spec,run_spec), z$n==n, z$p==p, z$seed==1000L+s)
  z
}
pick_rank <- function(curve) {
  if (length(curve)!=RMAX || any(!is.finite(curve))) return(NA_integer_)
  which.min(curve)
}
write_case <- function(records,n,p,prefix="") {
  selections <- do.call(rbind,lapply(records,function(z) z$selected))
  selections <- as.data.frame(selections)
  counts <- vapply(selections,function(v) as.integer(table(factor(v,levels=seq_len(RMAX)))),integer(RMAX))
  colnames(counts) <- names(selections)
  with_counts <- rbind(as.data.frame(counts),selections)
  export_dir <- file.path(result_dir,"reported_eps1e20")
  dir.create(file.path(export_dir,"raw_results"),recursive=TRUE,showWarnings=FALSE)
  stem <- case_stem(n,p,prefix)
  targets <- c(file.path(export_dir,paste0(stem,"_with_counts.rds")),
               file.path(export_dir,paste0(stem,"_with_counts.xlsx")),
               file.path(export_dir,"raw_results",paste0(stem,"_raw_results.rds")))
  # Export only once. If the replicate count changes, choose a new result_dir.
  if (all(file.exists(targets))) {
    if (!identical(readRDS(targets[1]),with_counts))
      stop("Existing export differs; choose a new result_dir: ",targets[1])
    return(invisible(targets))
  }
  if (any(file.exists(targets))) stop("Partial export retained; choose a new result_dir: ",stem)
  save_new(with_counts,targets[1])
  openxlsx::write.xlsx(with_counts,targets[2],overwrite=FALSE)
  save_new(records,targets[3])
  invisible(targets)
}

cv_fields <- function(prefix="") unlist(lapply(KS,function(k)
  paste0("facnum_",prefix,c("TCV","TPCV","TICV","TIPCV"),k,"K")),use.names=FALSE)
IC_NAMES <- c("IC1","IC2","zero","AIC1","AIC2","AIC3","BIC1","BIC2","BIC3",
              "IIC1","IIC2","IAIC1","IAIC2","IAIC3","IBIC1","IBIC2","IBIC3")
MAIN_FIELDS <- c("facnum_ratio","facnum_bcv","facnum_bcv2","facnum_sr",
                 paste0("facnum_",IC_NAMES),cv_fields(),cv_fields("NBknown_"),cv_fields("NBest_"))
