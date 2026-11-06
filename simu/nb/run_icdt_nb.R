
source("nb_settings.R")   # run from this folder
source(file.path(work_dir,"250113function_tcv.R"))   # generate_nb_data()
wrapper <- normalizePath(file.path(work_dir,"../../ICDT/icdt_wrapper.R"),mustWork=TRUE)
source(wrapper)

icdt_one <- function(s,n,p) {
  path <- rep_path(n,p,s,"icdt")
  if(file.exists(path)) {
    z <- check_record(readRDS(path),n,p,s)
    stopifnot(z$M==ICDT_M); return(z)
  }
  x <- generate_nb_data(n,p,s)$x
  error <- ""; warnings <- character(); started <- Sys.time()
  rr <- withCallingHandlers(tryCatch(
    run_icdt(x,qmax=RMAX,reptimes=ICDT_M,maxit=MAXITER,quiet=TRUE,restore_rng=FALSE),
    error=function(e) { error <<- conditionMessage(e); NULL }),
    warning=function(w) { warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning") })
  selected <- c(ICDT=NA_integer_,VSDT=NA_integer_,IC4=NA_integer_,IC5=NA_integer_)
  if(!is.null(rr)) {
    selected["ICDT"] <- pick_rank(rr$values$ICDT_val)
    selected["IC4"] <- pick_rank(rr$values$IC4_val)
    selected["IC5"] <- pick_rank(rr$values$IC5_val)
    if(all(is.finite(rr$values$VSDT_val_mat))) selected["VSDT"] <- rr$VSDT
  }
  out <- list(spec=run_spec,n=n,p=p,s=s,seed=1000L+s,M=ICDT_M,selected=selected,
              error=error,warnings=as.data.frame(table(warnings)),result=rr,
              elapsed_sec=as.numeric(difftime(Sys.time(),started,units="secs")),
              data_source="generate_nb_data(n,p,s), data seed 1000+s")
  save_new(out,path)
  cat("ICDT",n,p,"seed",1000L+s,"selected",selected["ICDT"],"\n")
  out
}

run_icdt_nb <- function() {
  cl <- NULL
  n_workers <- min(WORKERS,S)
  if(n_workers>1L) {
    log_path <- file.path(result_dir,sprintf("icdt_workers_%s_pid%d.log",
                                             format(Sys.time(),"%Y%m%d-%H%M%S"),Sys.getpid()))
    if(file.exists(log_path)) stop("Existing log retained: ",log_path)
    cl <- parallel::makePSOCKcluster(n_workers,outfile=log_path)
    on.exit(parallel::stopCluster(cl),add=TRUE)
    parallel::clusterCall(cl,function(wd,wr) {
      source(file.path(wd,"nb_settings.R"),local=globalenv())
      source(file.path(wd,"250113function_tcv.R"),local=globalenv())
      source(wr,local=globalenv()); TRUE
    },work_dir,wrapper)
    parallel::clusterExport(cl,"icdt_one",envir=environment())
  }
  for(n in N_VALUES) for(p in P_VALUES) {
    if(is.null(cl)) records <- lapply(seq_len(S),icdt_one,n=n,p=p)
    else records <- parallel::parLapplyLB(cl,seq_len(S),icdt_one,n=n,p=p,chunk.size=1L)
    write_case(records,n,p,"icdt_")
  }
}
if(sys.nframe()==0L) run_icdt_nb()
