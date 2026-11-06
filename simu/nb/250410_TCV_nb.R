
source("nb_settings.R")   # run from this folder
source(file.path(work_dir,"250113function_tcv.R"))

information_criteria <- function(obj,n,p,r) {
  # Same formulas and ordering as the original example.
  c(IC1=2*obj+r*((n+p)/(n*p))*log(n*p/(n+p)),
    IC2=2*obj+r/min(sqrt(n),sqrt(p))^2*log(min(sqrt(n),sqrt(p))^2),
    zero=2*obj,AIC1=2*obj+r*2/n,AIC2=2*obj+r*2/p,
    AIC3=2*obj+r*2*(n+p-r)/(n*p),BIC1=2*obj+r*log(n)/n,
    BIC2=2*obj+r*log(p)/p,BIC3=2*obj+r*(n+p-r)*log(n*p)/(n*p),
    IIC1=2*log(obj)+r*((n+p)/(n*p))*log(n*p/(n+p)),
    IIC2=2*log(obj)+r/min(sqrt(n),sqrt(p))^2*log(min(sqrt(n),sqrt(p))^2),
    IAIC1=2*log(obj)+r*2/n,IAIC2=2*log(obj)+r*2/p,
    IAIC3=log(2*obj)+r*2*(n+p-r)/(n*p),IBIC1=log(2*obj)+r*log(n)/n,
    IBIC2=log(2*obj)+r*log(p)/p,IBIC3=log(2*obj)+r*(n+p-r)*log(n*p)/(n*p))
}

run_one <- function(s,n,p) {
  path <- rep_path(n,p,s)
  if (file.exists(path)) return(check_record(readRDS(path),n,p,s))
  started <- Sys.time()
  cat(format(started),"START",n,p,"seed",1000L+s,"\n")
  warnings <- character(); errors <- list()
  z <- list(spec=run_spec,n=n,p=p,s=s,seed=1000L+s,
            selected=setNames(rep(NA_integer_,length(MAIN_FIELDS)),MAIN_FIELDS))
  attempt <- function(label,expr,otherwise=NA_real_) tryCatch(expr,error=function(e) {
    errors[[label]] <<- conditionMessage(e); otherwise
  })
  z <- withCallingHandlers({
    z <- c(z,generate_nb_data(n,p,s))
    x <- z$x
    # Baselines retain the preprocessing and settings of the previous examples.
    set.seed(1000L+s)
    z$selected["facnum_ratio"] <- attempt("Ratio",suppressMessages(
      chooseFacNumber_ratio(list(x),types="poisson",q_set=1:RMAX,select_method="SVR",verbose=FALSE)))
    z$selected["facnum_bcv"] <- attempt("BCV",{
      rr <- bcv::cv.svd.gabriel(sweep(log(x+.5),1,rowMeans(log(x+.5)),"-"),krow=2,kcol=2,maxrank=RMAX)
      as.numeric(names(which.min(colMeans(rr$msep)[-1])))
    })
    z$selected["facnum_bcv2"] <- attempt("CV",{
      rr <- bcv::cv.svd.wold(sweep(log(x+.5),1,rowMeans(log(x+.5)),"-"),k=5,maxrank=RMAX,tol=1e-3,maxiter=20)
      as.numeric(names(which.min(colMeans(rr$msep)[-1])))
    })
    z$selected["facnum_sr"] <- attempt("SR",{
      Rho <- cor(sweep(log(x+.5),1,rowMeans(log(x+.5)),"-"),method="spearman")
      SR(pracma::eig(Rho),Kmax=RMAX,p,n)
    })
    # Full-data fits also supply a data-driven pilot mean for size estimation.
    # No true rank or true mean is supplied to the pilot selector or size estimator.
    set.seed(2000L+s)
    full <- lapply(seq_len(RMAX),function(r) attempt(paste0("full_r",r),
      fit_poisson_nb_example(x,r,stage="full"),otherwise=NULL))
    objective <- vapply(full,function(f) if(is.null(f)) NA_real_ else f$nll,numeric(1))
    z$result_matrix <- t(vapply(seq_len(RMAX),function(r)
      information_criteria(objective[r],n,p,r),numeric(length(IC_NAMES))))
    colnames(z$result_matrix) <- paste0("log_",IC_NAMES)
    z$original_objective <- vapply(full,function(f) if(is.null(f)) NA_real_ else f$original_obj,numeric(1))
    z$full_fit_method <- vapply(full,function(f) if(is.null(f)) "failed" else f$method,character(1))
    z$full_am_error <- vapply(full,function(f) if(is.null(f)) "" else f$am_error,character(1))
    for (nm in IC_NAMES) z$selected[paste0("facnum_",nm)] <- pick_rank(z$result_matrix[,paste0("log_",nm)])
    pilot_rank <- pick_rank(z$result_matrix[,"log_IC2"])
    z$pilot_rank <- pilot_rank
    z$dispersion <- if(is.na(pilot_rank)) NULL else
      estimate_nb_size(x,full[[pilot_rank]]$mean)
    if (is.na(pilot_rank)) errors$pilot <- "Incomplete full-data JIC2 curve; estimated-size thinning unavailable."
    z$thinning <- list()
    for (k in KS) {
      pen <- seq_len(RMAX)*((n+p)/(n*p))*log((n*p)/(n+p))
      for (scheme in c("Poisson","NBknown","NBest")) {
        prefix <- if(scheme=="Poisson") "" else paste0(scheme,"_")
        sizes <- switch(scheme,Poisson=NULL,NBknown=z$true_size,
                        NBest=if(is.null(z$dispersion)) NULL else z$dispersion$size)
        if (scheme=="NBest" && is.null(z$dispersion)) next
        # Fixed streams make each scheme reproducible independently of loop order.
        set.seed(100000L+100L*s+k)
        rr <- attempt(paste0(scheme,"_K",k),multiDT(x,K=k,rmax=RMAX,overdisps=sizes),otherwise=NULL)
        if (is.null(rr)) next
        z$thinning[[paste0(scheme,"_K",k)]] <- rr
        curves <- list(TCV=rr$TCV,TPCV=rr$TCV+pen,TICV=rr$TICV,TIPCV=rr$TICV+pen)
        for (nm in names(curves)) {
          z[[paste0("log_",prefix,nm,k,"K")]] <- curves[[nm]]
          z$selected[paste0("facnum_",prefix,nm,k,"K")] <- pick_rank(curves[[nm]])
        }
      }
    }
    z
  },warning=function(w) {
    warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")
  })
  z$warnings <- as.data.frame(table(warnings),stringsAsFactors=FALSE)
  z$errors <- errors
  z$started <- started; z$finished <- Sys.time()
  z$elapsed_sec <- as.numeric(difftime(z$finished,started,units="secs"))
  z$status <- if(anyNA(z$selected)) "partial_failure" else "completed"
  save_new(z,path)
  cat(format(Sys.time()),"DONE",n,p,"seed",1000L+s,z$status,"seconds",round(z$elapsed_sec,1),"\n")
  z
}

run_main <- function() {
  cl <- NULL
  if(WORKERS>1L) {
    log_path <- file.path(result_dir,paste0("workers_",format(Sys.time(),"%Y%m%d-%H%M%S"),".log"))
    if(file.exists(log_path)) stop("Existing log retained: ",log_path)
    cl <- parallel::makePSOCKcluster(WORKERS,outfile=log_path)
    on.exit(parallel::stopCluster(cl),add=TRUE)
    doParallel::registerDoParallel(cl)
    parallel::clusterExport(cl,c("run_one","information_criteria"),envir=environment())
    parallel::clusterCall(cl,function(wd) {
      source(file.path(wd,"nb_settings.R"),local=globalenv())
      source(file.path(wd,"250113function_tcv.R"),local=globalenv())
      init_nb_cpp(); TRUE
    },work_dir)
  } else init_nb_cpp()
  for(n in N_VALUES) for(p in P_VALUES) {
    if(WORKERS>1L) {
      records <- foreach::foreach(s=seq_len(S)) %dopar% run_one(s,n,p)
    } else records <- lapply(seq_len(S),run_one,n=n,p=p)
    write_case(records,n,p)
  }
  session_path <- file.path(result_dir,"sessionInfo_main.txt")
  if(!file.exists(session_path)) capture.output(sessionInfo(),file=session_path)
}
if(sys.nframe()==0L) run_main()
