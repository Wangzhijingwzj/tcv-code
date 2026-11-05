# Counterpart of code/simu/simu-over/tcv_addmethods250506.R.
# Read the saved main-data matrix, so all methods use exactly the same data.
source("nb_settings.R")   # run from this folder
source(file.path(work_dir,"250113function_tcv.R"))

run_additional <- function() {
  for(n in N_VALUES) for(p in P_VALUES) {
    records <- lapply(seq_len(S),function(s) {
      path <- rep_path(n,p,s,"additional")
      if(file.exists(path)) return(check_record(readRDS(path),n,p,s))
      z <- check_record(readRDS(rep_path(n,p,s)),n,p,s)
      x <- z$x; errors <- list(); warnings <- character()
      attempt <- function(nm,expr) tryCatch(expr,error=function(e) {
        errors[[nm]] <<- conditionMessage(e); NA_real_
      })
      selected <- withCallingHandlers({
        log_centered_x <- sweep(log(x+.5),1,rowMeans(log(x+.5)),"-")
        eigvals_S <- attempt("cov_eigen",pracma::eig(cov(log_centered_x)))
        eigvals_R <- attempt("cor_eigen",pracma::eig(cor(log_centered_x,method="pearson")))
        eigvals_MKen <- attempt("Kendall_eigen",pracma::eig(SK(t(log_centered_x))))
        c(facnum_ne=attempt("SE_NE",NE(eigvals_S,p,n,maxfactor=RMAX)),
          facnum_ed=attempt("ED",ED(eigvals_S,Kmax=RMAX,4)),
          facnum_gr=attempt("GR",GR(eigvals_S,Kmax=RMAX,p,n)),
          facnum_er=attempt("ER",ER(eigvals_MKen,Kmax=RMAX,n,p)),
          facnum_act=attempt("ACT",ACT(eigvals_R,1+sqrt(p/n),Kmax=RMAX,n,p)))
      },warning=function(w) {
        warnings <<- c(warnings,conditionMessage(w)); invokeRestart("muffleWarning")
      })
      out <- list(spec=run_spec,n=n,p=p,s=s,seed=1000L+s,selected=selected,
                  errors=errors,warnings=as.data.frame(table(warnings)),
                  data_source=rep_path(n,p,s),
                  note="Original LADEL code is a constant 1 placeholder; it is not reported as a computed method.")
      save_new(out,path); out
    })
    write_case(records,n,p,"add_methods_")
  }
}
if(sys.nframe()==0L) run_additional()
