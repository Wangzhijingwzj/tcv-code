
# Based on code/simu/simu-over/250113function_tcv.R.
# Original helper functions below are retained verbatim. multiDT is extended
# only for NB allocation, the original +1e-20 validation loss, and fit diagnostics.
suppressPackageStartupMessages({
  library(foreach); library(doParallel); library(irlba); library(doSNOW)
  library(openxlsx); library(countsplit); library(bcv); library(tidyr)
  library(dplyr); library(invgamma); library(mvtnorm); library(pracma)
  library(pcaPP); library(Rfast)
})

init_nb_cpp <- function() {
  Rcpp::sourceCpp(file.path(work_dir,'confirm_jmle_omp_poisson_with_intercept_missing.cpp'),
    cacheDir=file.path(result_dir,'rcpp-cache',paste0('worker-',Sys.getpid())),
    cleanupCacheDir=FALSE,env=globalenv())
}

fit_poisson_nb_example <- function(X,r,stage=c('cv','full'),theta0=NULL,A0=NULL) {
  stage <- match.arg(stage); n <- nrow(X); p <- ncol(X)
  am_error <- ''
  fit <- tryCatch(suppressMessages(GFM::gfm(list(X),types='poisson',q=r,
     algorithm='AM',maxIter=MAXITER,verbose=FALSE)),error=function(e) {
       am_error <<- conditionMessage(e); NULL
     })
  if(!is.null(fit)) {
    mu <- exp(cbind(fit$hH,1)%*%t(cbind(fit$hB,fit$hmu)))
    if(all(is.finite(mu)) && all(mu>0)) {
      return(list(mean=mu,nll=mean(-dpois(X,mu,log=TRUE)),
                  original_obj=fit$obj,method='AM',am_error=am_error))
    }
    am_error <- 'AM returned nonfinite or nonpositive fitted means.'
  }
  if(is.null(theta0)) theta0 <- cbind(1,matrix(runif(n*r,0,1),n,r))
  if(is.null(A0)) A0 <- matrix(runif(p*(r+1),-1,1),p,r+1)
  fit <- confirm_CJMLE_poisson_cpp(X,matrix(TRUE,n,p),theta0,A0,
       matrix(TRUE,p,r+1),C=if(stage=='cv') Inf else 100,
       tol=if(stage=='cv') .01/(n*p) else .01)
  mu <- exp(fit$theta%*%t(fit$A))
  if(any(!is.finite(mu)) || any(mu<=0)) stop('CJMLE returned invalid fitted means.')
  # Normalize the objective for both solvers: GFM truncates probabilities,
  # whereas the original C++ obj omits the log-factorial term.
  list(mean=mu,nll=mean(-dpois(X,mu,log=TRUE)),original_obj=fit$obj,
       method='CJMLE',am_error=am_error)
}

multiDT <- function(x,K=2,rmax=8,overdisps=NULL) {
  n <- nrow(x); p <- ncol(x)
  xsplit <- countsplit::countsplit(x,folds=K,overdisps=overdisps)
  stopifnot(all(Reduce('+',lapply(xsplit,as.matrix))==x))
  obj_mat <- matrix(NA_real_,rmax,K)
  fit_method <- fit_error <- am_error <- matrix('',rmax,K)
  theta0 <- cbind(rep(1,n),matrix(runif(n*rmax,0,1),n,rmax))
  A0 <- matrix(runif(p*(rmax+1),-1,1),p,rmax+1)
  for(r in seq_len(rmax)) for(k in seq_len(K)) {
    XA <- x-as.matrix(xsplit[[k]]); XB <- x-XA
    fit <- tryCatch(fit_poisson_nb_example(XA,r,stage='cv',
      theta0=theta0[,seq_len(r+1),drop=FALSE],A0=A0[,seq_len(r+1),drop=FALSE]),
      error=function(e) { fit_error[r,k] <<- conditionMessage(e); NULL })
    if(is.null(fit)) { fit_method[r,k] <- 'failed'; next }
    lamB <- fit$mean/(K-1)
    obj_mat[r,k] <- mean(-log(dpois(XB,lamB)+1e-20))
    fit_method[r,k] <- fit$method; am_error[r,k] <- fit$am_error
  }
  cv <- rowSums(obj_mat)
  list(TCV=cv,TICV=log(cv),fold_loss=obj_mat,
       loss="Poisson_NLL_plus_1e-20",
       legacy_TCV=rowSums(obj_mat),legacy_fold_loss=obj_mat,
       fit_method=fit_method,fit_error=fit_error,am_error=am_error,
       size=overdisps,K=K)
}

# The same Gaussian factors, uniform loadings/intercepts, and three
# identifiability calls used in Example 3. Only the observation law changes.
generate_nb_data <- function(n,p,s) {
  set.seed(s+1000L)
  q <- TRUE_Q
  FF <- matrix(rnorm(n*q),nrow=n,ncol=q)
  BB <- matrix(runif(p*q,min=-1,max=1),nrow=p,ncol=q)
  a <- runif(p,min=0,max=1)
  FF0 <- add_identifiability(FF,BB,a)$H
  BB0 <- add_identifiability(FF,BB,a)$B
  alpha <- add_identifiability(FF,BB,a)$mu
  lambda <- exp(FF0%*%t(BB0)+matrix(1,n,1)%*%t(as.matrix(alpha)))
  true_size <- rep(NB_SIZE,p)
  x <- matrix(NA_real_,n,p)
  for(i in seq_len(n)) for(j in seq_len(p))
    x[i,j] <- rnbinom(1,size=true_size[j],mu=lambda[i,j])
  list(q=q,FF=FF,BB=BB,a=a,FF0=FF0,BB0=BB0,alpha=alpha,
       lambda=lambda,true_size=true_size,x=x)
}

estimate_nb_size <- function(x,mu) {
  # Feasible plug-in estimator: x and pilot fitted means only, no true values.
  bounds <- log(c(1e-4,1e6)); p <- ncol(x)
  ans <- lapply(seq_len(p),function(j) {
    ll <- function(log_size) sum(dnbinom(x[,j],mu=mu[,j],size=exp(log_size),log=TRUE))
    opt <- optimize(ll,interval=bounds,maximum=TRUE,tol=1e-7)
    candidates <- c(bounds,opt$maximum)
    vals <- vapply(candidates,ll,numeric(1))
    best <- which.max(vals)
    c(size=exp(candidates[best]),loglik=vals[best],boundary=as.integer(best<=2L))
  })
  ans <- as.data.frame(do.call(rbind,ans))
  list(size=ans$size,table=ans,bounds=exp(bounds))
}

add_identifiability <- function(H, B, mu){
  # Load the irlba library
  #library(irlba)
  
  # Perform SVD decomposition with rank k = 10
  
  mu <- mu + B %*% colMeans(H)
  q <- ncol(H); n <- nrow(H)
  svdHB <- irlba((H- matrix(colMeans(H), n, q, byrow = TRUE)) %*% t(B), nv= q)
  #signB1 <- sign(svdHB$v[1,])
  signB1 = sign(apply(svdHB$v, 2, function(x) x[which.max(abs(x))]))
  H <- sqrt(n) * svdHB$u %*% Diag(signB1)
  
  B <- svdHB$v %*% Diag(svdHB$d[1:q]*signB1) / sqrt(n)
  
  return(list(H=H, B=B, mu=mu))
}




chooseFacNumber_ratio <- function(XList, types, q_set = 2: 10, select_method = c("SVR", "IC"),offset=FALSE, dc_eps=1e-4, maxIter=30,
                                  verbose = TRUE, parallelList=NULL){
  
  
  q_max <- max(q_set)
  gfm2 <- GFM::gfm(XList, types, q=q_max, algorithm = 'AM', offset=offset, dc_eps=dc_eps, maxIter=maxIter,verbose = verbose)
  
  
  svalues <- svd(gfm2$hB)$d
  svalues_use <- svalues[svalues>1e-2]
  q_max <- length(svalues_use)
  q <- which.max(svalues[-q_max] / svalues[-1])
  
  message('The singular value ratio (SVR) method estimates the factor number q  as ', q, '. \n')
  
  return(q)
  
}




SR <- function(eigvals, Kmax, dim, size) {
  m_prime_value <- rep(0, Kmax + 1)
  for (j in 1:(Kmax + 1)) {
    cj <- (dim - j) / size
    m_prime_value[j] <- mean(na.omit(1 / ((eigvals[-(1:j)] - eigvals[j])^2+dim^{-4/3}) ) )
  }
  m_prime_ratio <- m_prime_value[1:Kmax + 1] / m_prime_value[1:Kmax]
  Khat_SR <- which.max(m_prime_ratio)
  return(Khat_SR)
}


# ------------------------------------------------------------------------------------
# Function: Estimate population spiked eigenvalues from sample eigenvalues
# Reference: Bai and Ding (2012, RMTA), Fan et al. (2020, JASA)
# Input:
#         1. eigvals: sample eigenvalues of Pearson's correlation matrix,
#         2. Kmax: predetermined upper bound of the true number of spiked eigenvalues
# Output:
#         1. eigvals_C: Kmax estimated spiked eigenvalues
# ------------------------------------------------------------------------------------
BaiDing12_RMTA <- function(eigvals, Kmax) {
  mz <- rep(0, Kmax)
  for (kk in 1:Kmax) {
    qu <- 3 / 4
    eigvals_1 <- eigvals[-seq(1, kk, 1)]
    z0 <- qu * eigvals[kk] + (1 - qu) * eigvals[kk + 1]
    ttt <- c(eigvals_1 - eigvals[kk], z0 - eigvals[kk])
    y0 <- (length(eigvals_1)) / (n - 1)
    mz[kk] <- -(1 - y0) / eigvals[kk] + y0 * (sum(na.omit(1 / ttt)) / (p - kk))
  }
  eigvals_C <- -1 / mz
  return(eigvals_C)
}
# ------------------------------------------------------------------------
# Function: ACT method
# Reference: Fan et al. (2020, JASA)
# Input:
#         1. eigvals: sample eigenvalues of Pearson's correlation matrix
#         2. threshold: threshold for adjusted eigenvalues
#         3. Kmax: predetermined upper bound of the true number of factors
# Output:
#         1. Khat_ACT: estimated number of common factors
# ------------------------------------------------------------------------
ACT <- function(eigvals, threshold, Kmax, n ,p) {
  mz <- rep(0, Kmax)
  for (kk in 1:Kmax) {
    qu <- 3 / 4
    eigvals_1 <- eigvals[-seq(1, kk, 1)]
    z0 <- qu * eigvals[kk] + (1 - qu) * eigvals[kk + 1]
    ttt <- c(eigvals_1 - eigvals[kk], z0 - eigvals[kk])
    y0 <- (length(eigvals_1)) / (n - 1)
    mz[kk] <- -(1 - y0) / eigvals[kk] + y0 * (sum(na.omit(1 / ttt)) / (p - kk))
  }
  eigvals_C <- -1 / mz
  #eigvals_C <- BaiDing12_RMTA(eigvals, Kmax)
  temp <- eigvals_C - threshold
  temp1 <- seq(1, Kmax, 1)
  temp2 <- cbind(temp1, temp)
  Khat_ACT <- max((temp2[, 1][temp2[, 2] > 0]), 0)
  return(Khat_ACT)
}


getCurrentPath <- function() {
  this_file <- commandArgs() %>%
    tibble::enframe(name = NULL) %>%
    tidyr::separate(col = value, into = c("key", "value"), sep = "=", fill = "right") %>%
    dplyr::filter(key == "--file") %>%
    dplyr::pull(value)
  # If script is not run from command line, get path from RStudio
  if (length(this_file) == 0) {
    this_file <- rstudioapi::getSourceEditorContext()$path
  }
  return(dirname(this_file))
}


library(invgamma)
library(mvtnorm)
library(pracma)
library(pcaPP)
library(Rfast)


# -------------------------------------------------------------------
# Function: Calculate Multivariate Kendall's tau matrix (MKen)
# Input:
#         1. X: p x n, data matrix (p = dimension, n = sample size)
# Output:
#         1. TK: p x p, MKen of X
# -------------------------------------------------------------------
SK <- function(X) {
  p <- nrow(X)
  n <- ncol(X)
  TK <- matrix(0, p, p)
  for (i in 1:(n - 2)) {
    TT <- X[, i] - X[, (i + 1):n]
    V <- sqrt(colsums(TT^2))
    TT_spatial_sign <- sweep(TT, 2, V, FUN = "/")
    TT <- tcrossprod(TT_spatial_sign)
    TK <- TK + TT
  }
  TT <- X[, n - 1] - X[, n]
  TK <- TK + TT %*% t(TT) / sum(TT^2)
  TK <- 2 / (n * (n - 1)) * TK
  return(TK)
}



# ---------------------------------------------------------------------------
# Function: ER method
# Reference: Lam and Yao (2012, AoS), Ahn and Horenstein (2013, Econometrica)
# Input:
#         1. eigvals: sample eigenvalues of some kinds of matrices, e.g.,
#                     sample covariance matrix
#         2. Kmax: predetermined upper bound of the true number of factors
# Output:
#         1. Khat_ER: estimated number of common factors
# ---------------------------------------------------------------------------
ER <- function(eigvals, Kmax,n ,p) {
  eigvals_lag1 <- c(eigvals[-1], 0)
  ER <- (eigvals_lag1 / eigvals)[-n]
  Khat_ER <- which.min(ER[1:Kmax])
  return(Khat_ER)
}


# ---------------------------------------------------------------------------
# Function: GR method
# Reference: Ahn and Horenstein (2013, Econometrica)
# Input:
#         1. eigvals: sample eigenvalues of some kinds of matrices, e.g.,
#                     sample covariance matrix
#         2. Kmax: predetermined upper bound of the true number of factors
# Output:
#         1. Khat_GR: estimated number of common factors
# ---------------------------------------------------------------------------
GR <- function(eigvals, Kmax, dim, size) {
  hm1 <- hm2 <- hm3 <- matrix(rep(0, dim * Kmax), dim, Kmax)
  for (i in 1:Kmax) {
    for (j in i:dim) {
      hm1[j, i] <- 1
    }
  }
  for (i in 1:Kmax) {
    for (j in (i + 1):dim) {
      hm2[j, i] <- 1
    }
  }
  for (i in 1:Kmax) {
    for (j in (i + 2):dim) {
      hm3[j, i] <- 1
    }
  }
  gr1 <- eigvals %*% hm1
  gr2 <- eigvals %*% hm2
  gr3 <- eigvals %*% hm3
  GR <- log(gr1 / gr2) / log(gr2 / gr3)
  Khat_GR <- which.max(GR)
  return(Khat_GR)
}



# ------------------------------------------------------------------------
# Function: ED method
# Reference: Onatski (2010, RES)
# Input:
#         1. eigvals: sample eigenvalues of sample covariance matrix
#         2. Kmax: predetermined upper bound of the true number of factors
#         3. iter: number of iteration
# Output:
#         1. Khat_ED: estimated number of common factors
# ------------------------------------------------------------------------
ED <- function(eigvals, Kmax, iter) {
  iter <- 4
  ED_seq <- eigvals[1:Kmax] - eigvals[2:(Kmax + 1)]
  noj <- Kmax + 1
  for (j in 1:iter) {
    y <- eigvals[noj:(noj + iter)]
    x <- (noj + seq(-1, (iter - 1), 1))^(2 / 3)
    beta <- sum((x - mean(x)) * (y - mean(y))) / sum((x - mean(x))^2)
    delta <- 2 * abs(beta)
    # min factor number is 1
    # Khat_ED <- ifelse(max(ED_seq - delta) < 0, 0, max(which(ED_seq >= delta)))
    Khat_ED <- ifelse(max(ED_seq - delta) < 0, 1, max(which(ED_seq >= delta)))
    noj <- Khat_ED + 1
  }
  return(Khat_ED)
}


# ------------------------------------------------------------------------
# Function: NE method
# Reference: Nadakuditi and Edelman (2008, IEEE Transactions on Signal Processing)
# Input:
#         1. eigvals: sample eigenvalues of Pearson's correlation matrix
#         2. Kmax: predetermined upper bound of the true number of factors
#         3. dim: dimension
#         4. size: sample size
# Output:
#         1. Khat_SR: estimated number of common factors
# ------------------------------------------------------------------------
NE <- function(eigvals, dim, size, maxfactor=min(p, n)) {
  p <- dim
  n <- size
  t <- rep(0, maxfactor)
  ## min number of factor is one
  for (i in 1:(maxfactor)) {
    temp <- p * ((p - i) * (sum(eigvals[(i + 1):p]^2) / sum(eigvals[(i + 1):p])^2) - 1 - p / n) - p / n
    t[i] <- temp^2 * (n / p)^2 / 4 + 2 * (i + 1)
  }
  Khat_NE <- which.min(t)
  # return(t)
  return(Khat_NE)
}

