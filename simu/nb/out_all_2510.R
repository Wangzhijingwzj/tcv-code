# Summary under the original example's -log(dpois()+1e-20) validation loss.
# Older checkpoints are rescored from stored fold losses, without refitting.
source("nb_settings.R")   # run from this folder

summarize_nb <- function() {
  rows <- list(); dispersion <- list(); fits <- list(); issues <- list(); cases <- list()
  for(n in N_VALUES) for(p in P_VALUES) for(s in seq_len(S)) {
    for(stage in c("main","additional","icdt")) {
      path <- rep_path(n,p,s,stage)
      if(!file.exists(path)) next
      z <- check_record(readRDS(path),n,p,s)
      key <- paste(n,p,stage,sep="_")
      cases[[key]] <- c(cases[[key]],list(z))
      rows[[length(rows)+1L]] <- data.frame(n=n,p=p,seed=1000L+s,stage=stage,
        method=names(z$selected),selected=as.numeric(z$selected),source=path)
      if(stage=="main") {
        if(!is.null(z$dispersion)) dispersion[[length(dispersion)+1L]] <-
          data.frame(n=n,p=p,seed=z$seed,pilot_rank=z$pilot_rank,taxon=seq_len(p),
                     true_size=z$true_size,z$dispersion$table)
        for(nm in names(z$thinning)) {
          rr <- z$thinning[[nm]]
          grid <- expand.grid(rank=seq_len(RMAX),fold=seq_len(rr$K))
          fits[[length(fits)+1L]] <- data.frame(n=n,p=p,seed=z$seed,scheme=nm,grid,
            method=as.vector(rr$fit_method),loss=as.vector(rr$fold_loss),
            am_error=as.vector(rr$am_error),error=as.vector(rr$fit_error))
        }
      }
      if(!is.null(z$warnings) && nrow(z$warnings)) issues[[length(issues)+1L]] <-
        data.frame(n=n,p=p,seed=z$seed,stage=stage,type="warning",
                   message=as.character(z$warnings[[1]]),count=z$warnings[[2]])
      errs <- if(stage=="icdt") z$error else unlist(z$errors,use.names=TRUE)
      if(length(errs)) for(e in errs[nzchar(errs)]) issues[[length(issues)+1L]] <-
        data.frame(n=n,p=p,seed=z$seed,stage=stage,type="error",message=e,count=1)
    }
  }
  if(!length(rows)) stop("No saved replicates to summarize.")
  # Regenerate complete per-case exports in a new directory; old exports stay intact.
  for(key in names(cases)) if(length(cases[[key]])==S) {
    zz <- cases[[key]]; stage <- sub("^[0-9]+_[0-9]+_","",key)
    prefix <- switch(stage,main="",additional="add_methods_",icdt="icdt_")
    write_case(zz,zz[[1]]$n,zz[[1]]$p,prefix)
  }
  runs <- do.call(rbind,rows)
  summary <- do.call(rbind,lapply(split(runs,interaction(runs$n,runs$p,runs$method,drop=TRUE)),function(d) {
    v <- d$selected[is.finite(d$selected)&d$selected>=1&d$selected<=RMAX&d$selected==floor(d$selected)]
    co <- sum(v==TRUE_Q); ov <- sum(v>TRUE_Q); un <- sum(v<TRUE_Q)
    data.frame(n=d$n[1],p=d$p[1],method=d$method[1],requested=S,attempted=nrow(d),
      missing=S-nrow(d),valid=length(v),failed=nrow(d)-length(v),correct=co,over=ov,under=un,
      MO=if(ov) mean(v[v>TRUE_Q]-TRUE_Q) else 0,MU=if(un) mean(TRUE_Q-v[v<TRUE_Q]) else 0,
      correct_over_under=sprintf("%d (%d|%d)",co,ov,un),
      mean_over_under=sprintf("%.2f|%.2f",if(ov) mean(v[v>TRUE_Q]-TRUE_Q) else 0,
                             if(un) mean(TRUE_Q-v[v<TRUE_Q]) else 0))
  })); rownames(summary) <- NULL
  paper_names <- c("facnum_bcv","facnum_bcv2","facnum_ne","facnum_sr","facnum_ratio",
    "facnum_IC1","facnum_IC2","ICDT","VSDT",cv_fields(),cv_fields("NBknown_"),cv_fields("NBest_"),
    "facnum_ed","facnum_gr","facnum_er","facnum_act")
  paper_names <- paper_names[!grepl("_TICV|^facnum_TICV",paper_names)]
  comparison <- summary[summary$method %in% paper_names,]
  comparison <- comparison[order(comparison$n,comparison$p,match(comparison$method,paper_names)),]
  wide_names <- c("n","p",unlist(lapply(paper_names,function(m)
    paste0(m,c("_correct(over|under)","_over_mean|under_mean"))),use.names=FALSE))
  wide <- do.call(rbind,lapply(split(comparison,interaction(comparison$n,comparison$p,drop=TRUE)),function(d) {
    out <- setNames(as.list(rep(NA_character_,length(wide_names))),wide_names)
    out$n <- d$n[1]; out$p <- d$p[1]
    for(i in seq_len(nrow(d))) {
      out[[paste0(d$method[i],"_correct(over|under)")]] <- d$correct_over_under[i]
      out[[paste0(d$method[i],"_over_mean|under_mean")]] <- d$mean_over_under[i]
    }
    as.data.frame(out,check.names=FALSE)
  }))
  dest <- file.path(result_dir,paste0("output_eps1e20_",format(Sys.time(),"%Y%m%d-%H%M%S")))
  if(dir.exists(dest)) stop("Existing summary retained: ",dest)
  dir.create(dest)
  bind_or_empty <- function(x) if(length(x)) do.call(rbind,x) else data.frame()
  notes <- data.frame(note=c(
    sprintf("NB data with true size=%g and rank=%d; completed checkpoints may be reused.",NB_SIZE,TRUE_Q),
    "All TCV/TPCV/TIPCV comparisons use mean(-log(dpois(XB,lamB)+1e-20)), matching the original examples.",
    "Old checkpoints are rescored from stored +1e-20 fold losses. Data, fitted models and nuisance-size estimates are not recomputed.",
    "The full-data pilot/JIC implementation and other comparison methods retain their existing definitions; only TCV validation scoring changed.",
    "NBknown uses true size; NBest estimates taxon-specific sizes using its existing JIC2 pilot. Estimated size does not ensure exact independence.",
    "Exact-loss fields in historical raw records are archival diagnostics and are not used in this report.",
    "Completed per-case exports under the requested score are in reported_eps1e20/; original exports remain unchanged.",
    "ICDT uses the existing author implementation with M=3 and its original fallback.",
    "Failures and missing observations are reported separately. Finite results do not certify convergence."))
  tabs <- list(paper_comparison=comparison,paper_format=wide,all_methods=summary,replicates=runs,
    dispersion=bind_or_empty(dispersion),fit_diagnostics=bind_or_empty(fits),
    warnings_errors=bind_or_empty(issues),notes=notes)
  save_new(tabs,file.path(dest,"all_result2510_nb.rds"))
  openxlsx::write.xlsx(tabs,file.path(dest,"all_result2510_nb.xlsx"),overwrite=FALSE)
  write.csv(comparison,file.path(dest,"comparison.csv"),row.names=FALSE)
  print(comparison[,c("n","p","method","valid","failed","correct","over","under","MO","MU")],row.names=FALSE)
  cat("SUMMARY:",dest,"\n")
  invisible(dest)
}
if(sys.nframe()==0L) summarize_nb()
