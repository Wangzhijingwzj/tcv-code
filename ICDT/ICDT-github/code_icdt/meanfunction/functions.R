
################################################################################
################ functions needed for ICDT
################################################################################

# compute the mode
getModeAverage <- function(v) {
  tab <- table(v)
  modes <- names(tab[tab == max(tab)])
  if(length(modes) > 1) {
    modes <- as.numeric(modes)
    modeAverage <- mean(modes)
  } else {
    if(is.numeric(v)) modes <- as.numeric(modes)
    modeAverage <- modes
  }
  return(modeAverage)
}

# check the dataset to ensure that no row or column consists entirely of identical values
check_mat <- function(x) {
  for (i in 1:nrow(x)) {
    if (length(unique(x[i, ])) == 1) {
      cols <- sample(1:ncol(x), size = 5, replace = FALSE)
      x[i, cols] <- x[i, cols] + 2
    }
  }
  
  for (j in 1:ncol(x)) {
    if (length(unique(x[, j])) == 1) {
      rows <- sample(1:nrow(x), size = 5, replace = FALSE)
      x[rows, j] <- x[rows, j] + 2
    }
  }
  
  return(x)
}

# Impute missing values (NA) using column-wise means
replace_with_column_means <- function(H) {
  for (j in 1:ncol(H)) {
    column_mean <- mean(H[,j][!is.na(H[,j]) & !is.nan(H[,j]) & !is.infinite(H[,j])], na.rm = TRUE)
    H[,j][is.na(H[,j]) | is.nan(H[,j]) | is.infinite(H[,j])] <- column_mean
  }
  return(H)
}


# function for element splitting 
split_matrix <- function(n, p, K) {
  original_matrix <- matrix(1, nrow = n, ncol = p)
  split_matrices <- replicate(K, matrix(0, nrow = n, ncol = p), simplify = FALSE)
  for (i in 1:n) {
    for (j in 1:p) {
      counts <- rmultinom(1, size = original_matrix[i, j], prob = rep(1/K, K))
      for (k in 1:K) {
        split_matrices[[k]][i, j] <- counts[k]
      }
    }
  }
  split_matrices <- lapply(split_matrices, function(x) x == 1)
  return(split_matrices)
}

# function for computing the likelihood
JL <- function(resp, theta, A, nonmis_ind){
  temp = theta%*% t(A)
  temp1 = resp * temp - exp(temp)  
  lik = sum(temp1[nonmis_ind])
  list(lik = lik, M = temp)
}

# The following gfm functions and their dependencies are derived from the R package GFM 
# and have been reimplemented here for debugging purposes
# see ''Generalized Factor Model for Ultra-High Dimensional Correlated Variables with Mixed Types'' (JASA 2023) for more details


Diag <- function(vec){
  q <- length(vec)
  if(q > 1){
    y <- diag(vec)
  }else{
    y <- matrix(vec, 1,1)
  }
  return(y)
}

cor.mat <- function (p, rho, type = "toeplitz"){
  if(p == 1) return(matrix(1,1,1))
  mat <- diag(p)
  if (type == "toeplitz") {
    for (i in 2:p) {
      for (j in 1:i) {
        mat[i, j] <- mat[j, i] <- rho^(abs(i - j))
      }
    }
  }
  if (type == "identity") {
    mat[mat == 0] <- rho
  }
  return(mat)
}

gendata <- function (seed = 1, n = 300, p = 50,
                     type =  c('homonorm', 'heternorm', 'pois', 'bino', 'norm_pois', 'pois_bino', 'npb'),
                     q = 6, rho = 1, n_bin=1){
  # library(MASS)
  # Diag <- GFM:::Diag
  # cor.mat <- GFM:::cor.mat
  
  type <- match.arg(type)
  if (length(rho) == 1)
    rho <- c(rho, 1.5)
  factor_term <- rho[2]
  
  Z <- matrix(rnorm(p * q), p, q)
  B <- qr(Z)
  eigvs <- sqrt(sort(eigen(t(Z) %*% Z)$values, decreasing = T))
  B1 <- qr.Q(B) %*% Diag(sqrt(eigvs))
  B0 <- rho[1] * B1 %*% Diag(sign(B1[1, ]))
  mu0 <- 0.4 * rnorm(p)
  Bm0 <- cbind(mu0, B0)
  set.seed(seed)
  H <- mvrnorm(n, mu = rep(0, q), cor.mat(q, 0.5))
  svdH <- svd(cov(H))
  H0 <- scale(H, scale = F) %*% svdH$u %*% Diag(1/sqrt(svdH$d)) %*%
    svdH$v
  if (type == "homonorm") {
    X <- H0 %*% t(B0) + matrix(mu0, n, p, byrow = T) + mvrnorm(n,
                                                               rep(0, p), diag(p))
    group <- rep(1, p)
    XList <- list(X)
    types <- c("gaussian")
    
  }else if (type == "heternorm") {
    sigmas = 0.1 + 4 * runif(p)
    X <- H0 %*% t(B0) + matrix(mu0, n, p, byrow = T) + mvrnorm(n,
                                                               rep(0, p), diag(sigmas))
    group <- rep(1, p)
    
    XList <- list(X)
    types <- c("gaussian")
    
  }else if (type == "pois") {
    g1 <- 1:p
    B0[g1, ] <- B0[g1, ]/max(B0[g1, ]) * factor_term
    mu <- exp(H0 %*% t(B0) + matrix(mu0, n, p, byrow = T))
    X <- matrix(rpois(n * p, lambda = mu), n, p)
    group <- rep(1, p)
    XList <- list(X[,g1])
    types <- c("poisson")
  }else if(type == 'bino'){
    g1 <- 1:p
    mu <- 1/(1 + exp(-cbind(1, H0) %*% t(Bm0[g1, ])))
    X <- matrix(rbinom(prod(dim(mu)), n_bin, mu), n, p)
    group <- rep(1, p)
    
    XList <- list(X[,g1])
    types <- c("binomial")
  }else if (type == "norm_pois") {
    g1 <- 1:floor(p/2)
    g2 <- (floor(p/2) + 1):p
    Bm0[g2, -1] <- Bm0[g2, -1]/max(Bm0[g2, -1]) * factor_term
    mu1 <- cbind(1, H0) %*% t(Bm0[g1, ])
    mu2 <- exp(cbind(1, H0) %*% t(Bm0[g2, ]))
    X <- cbind(matrix(rnorm(prod(dim(mu1)), mu1, 1), n, floor(p/2)),
               matrix(rpois(prod(dim(mu2)), mu2), n, ncol(mu2)))
    group <- c(rep(1, length(g1)), rep(2, length(g2)))
    
    XList <- list(X[,g1], X[,g2])
    types <- c("gaussian", "poisson")
    
  }else if (type == "pois_bino") {
    g1 <- 1:floor(p/2)
    g2 <- (floor(p/2) + 1):p
    Bm0[g1, -1] <- Bm0[g1, -1]/max(Bm0[g1, -1]) * factor_term
    mu1 <- exp(cbind(1, H0) %*% t(Bm0[g1, ]))
    mu2 <- 1/(1 + exp(-cbind(1, H0) %*% t(Bm0[g2, ])))
    X <- cbind(matrix(rpois(prod(dim(mu1)), mu1), n, ncol(mu1)),
               matrix(rbinom(prod(dim(mu2)), n_bin, mu2), n, ncol(mu2)))
    group <- c(rep(1, length(g1)), rep(2, length(g2)))
    XList <- list(X[,g1], X[,g2])
    types <- c("poisson", 'binomial')
  }else if(type == 'npb'){
    g1 <- 1:floor(p/3)
    g2 <- (floor(p/3) + 1):floor(p*2/3)
    g3 <- (floor(2*p/3) + 1):p
    mu11 <- cbind(1, H0) %*% t(Bm0[g1, ])
    Bm0[g2, -1] <- Bm0[g2, -2]/max(Bm0[g2, -1]) * factor_term
    mu1 <- exp(cbind(1, H0) %*% t(Bm0[g2, ]))
    mu2 <- 1/(1 + exp(-cbind(1, H0) %*% t(Bm0[g3, ])))
    X <- cbind(matrix(rnorm(prod(dim(mu11)),mu11, 1), n, ncol(mu11)),
               matrix(rpois(prod(dim(mu1)), mu1), n, ncol(mu1)),
               matrix(rbinom(prod(dim(mu2)), n_bin, mu2), n, ncol(mu2)))
    group <- c(rep(1, length(g1)), rep(2, length(g2)), rep(3, length(g3)))
    XList <- list(X[,g1], X[,g2], X[,g3])
    types <- c("gaussian", "poisson", 'binomial')
  }
  
  return(list(X=X, XList = XList,  types= types, B0 = B0, H0 = H0, mu0 = mu0))
}




gfm <- function(XList, types, q=10, offset=FALSE, dc_eps=1e-4, maxIter=30,
                verbose = FALSE, algorithm=c("VEM", "AM")){
  
  
  if((!is.null(q)) && (q<1) ) stop("q must be NULL or other positive integer!")
  n <- nrow(XList[[1]]); p <- sum(sapply(XList, ncol))
  if(p <20) stop("ncol(X) must be at least no less than 20!")
  if(n <20) stop("nrow(X) must be at least no less than 20!")
  
  algorithm <- match.arg(algorithm)
  
  type_map <- 1:3;
  names(type_map) <- c("gaussian", "poisson", "binomial")
  typeID <- unname(type_map[types]) ## match IDs
  if(length(setdiff(types,names(type_map)))>1){
    stop("The component of string vector types must be contained in ('gaussian', 'poisson', 'binomial')!")
  }
  
  ### two-step estimation
  ### algorithm 1: use VEM to obtain initial values, then use one-step update to improve the estimate.
  ### algorithm 2: use alternate maximization to obtain initial values, then use one-step update to improve the estimate.
  if(algorithm== "VEM"){
    if(verbose){
      message('Starting the two-step method with varitional EM in the first step...\n')## use VEM to obtain initial values
    }
    
    tic <- proc.time()
    reslist <- vem.fit(XList, types, q=q, offset=offset, epsELBO=dc_eps, maxIter=maxIter, verbose=verbose)
    toc <- proc.time()
    if(verbose){
      message('Finish the two-step method\n')
    }
    
    gfm2 <- list()
    gfm2$hH <- reslist$H
    gfm2$hB <- reslist$B
    gfm2$hmu <- t(reslist$mu)
    gfm2$obj <- reslist$ELBO
    gfm2$history <- list(c=reslist$ELBO_seq, maxIter=maxIter, eplasedTime=toc-tic)
    
    try({
      X <- NULL
      group <- NULL
      for(i in seq_along(typeID)){
        X <- cbind(X, XList[[i]])
        group <- c(group, rep(i, ncol(XList[[i]])))
      }
      rm(XList)
      gfm2final <- gfm_eval_intercept_osfinal(X, gfm2$hH, gfm2$hB,gfm2$hmu, group, types)
      if(any(is.nan(gfm2final$hH)) || any(is.nan(gfm2final$hB))) stop("there is NaNs produced!")
      gfm2final$history <- gfm2$history
      gfm2 <- gfm2final
    }, silent = TRUE)
    
  }else if(algorithm== 'AM'){
    
    X <- NULL
    group <- NULL
    for(i in seq_along(typeID)){
      X <- cbind(X, XList[[i]])
      group <- c(group, rep(i, ncol(XList[[i]])))
    }
    rm(XList)
    omega <- 1/p
    if(verbose){
      message('Starting the two-step method with alternate maximization in the first step...\n')
    }
    
    dropout <- 0
    if(any(typeID==2) && length(typeID)>1) dropout <- which(typeID==2)
    gfm2 <- gfm_eval_intercept_init(X, group, types, q,
                                    dropout, dc_eps, maxIter,
                                    verbose)
    
    if(verbose){
      message('Finish the two-step method\n')
    }
    
    try({
      gfm2final <- gfm_eval_intercept_osfinal(X, gfm2$hH, gfm2$hB,gfm2$hmu, group, types)
      if(any(is.nan(gfm2final$hH)) || any(is.nan(gfm2final$hB))) stop("there is NaNs produced!")
      gfm2final$history <- gfm2$history
      gfm2 <- gfm2final
    }, silent = TRUE)
  }
  ## Add identifiable condition
  res_idents <- add_identifiability(gfm2$hH, gfm2$hB, gfm2$hmu)
  gfm2$hH <- res_idents$H
  gfm2$hB <- res_idents$B
  gfm2$hmu <- res_idents$mu
  gfm2$q <- q
  class(gfm2) <- 'gfm'
  return(gfm2)
}


singleIC <- function(X, group, type, q_set=1:10, dropout=0, dc_eps=1e-4, maxIter=10,output=FALSE){
  n <- nrow(X)
  q_num <- length(q_set)
  allhH <- list(); allhB <- list()
  for(r in 1:q_num){
    
    gfm1 <- gfm_eval_intercept_init(X, group, type, r,
                                    dropout, dc_eps, maxIter, output)
    
    
    gfm1 <-gfm_eval_intercept_osfinal(X, gfm1$hH, gfm1$hB,gfm1$hmu, group, type)
    
    allhH[[r]] <- cbind(1, gfm1$hH)
    allhB[[r]] <- cbind(gfm1$hmu, gfm1$hB)
  }
  Vr <- matrix(0, 2, q_num)
  for(r in 1:q_num){
    Vr[,r] <- ICriteria(X, allhB[[r]], allhH[[r]], r, group, type)
    
  }
  IC <- apply(Vr, 2, sum)
  q <- q_set[which.min(IC)]
  return(q)
}



ICriteria <- function(X, hB, hH, r, group, type, criteria){
  # hB <- allhB[[r]]; hH <- allhH[[r]]
  n <- nrow(X); p <- ncol(X)
  omega <- 1/p
  ind_set <- unique(group)
  ng <- length(ind_set)
  gcell <- list()
  for(j in 1:ng){
    gcell[[j]] <- which(group ==j)
  }
  
  c1 <- objfunc(hH, hB, X, omega, gcell, type)
  Vr <- switch(criteria,
               IC=c(log(c1+1e-7), r/min(sqrt(n), sqrt(p))^2*log(min(sqrt(n), sqrt(p))^2)),
               PC=c(c1, r/min(sqrt(n), sqrt(p))^2*log(min(sqrt(n), sqrt(p))^2)),
               IC1=c(c1, r*((n+p)/(n*p))*log(n*p/(n+p))),
               IC2=c(c1, r/min(sqrt(n), sqrt(p))^2*log(min(sqrt(n), sqrt(p))^2)),
               zero=c(c1, 0),
               AIC1=c(c1,r*2/n),
               AIC2=c(c1,r*2/p),
               AIC3=c(c1,r*2*(n+p-r)/(n*p)),
               BIC1=c(c1, r*log(n)/n),
               BIC2=c(c1, r*log(p)/p),
               BIC3=c(c1, r*(n+p-r)*log(n*p)/(n*p)),
               pen=c(c1, (r*n+(r+1)*p)/(n*p)*log((n*p)/(r*n+(r+1)*p))),
               # r * (n+p)/(n*p)*log(n*p/(n+p))
               pen1=c(c1, (r*n+(r+1)*p)/(n*p)*0.5*log((n*p)/(n+p))),
               pen2=c(c1, (r*n+(r+1)*p)/(n*p)*(0.5*log((n*p)/(n+p))-1)),
               pen3=c(c1, (r*n+(r+1)*p)/(n*p)*log((0.5*n*p)/(n+p))),
               pen4=c(c1, (r*n+(r+1)*p)/(n*p)*(log((0.5*n*p)/(n+p))-1)),
               pen5=c(c1, (r*n+r*p)/(n*p)*0.5*log((n*p)/(n+p))),
               pen6=c(c1, (r*n+r*p)/(n*p)*(0.5*log((n*p)/(n+p))-1)),
               pen7=c(c1, (r*n+r*p)/(n*p)*log((0.5*n*p)/(n+p))),
               pen8=c(c1, (r*n+r*p)/(n*p)*(log((0.5*n*p)/(n+p))-1)),
               cov1=c(c1, log((n*p)/(n+p))),
               cov2=c(c1 , log(min(sqrt(n), sqrt(p))^2))
               
  )
  return(Vr)
}


paraIC <- function(r, XX, group, type, cri,
                   dropout=0, eps2=1e-4, maxIter=10, output=FALSE, fast_version=TRUE){
  
  gfm1 <- gfm_eval_intercept_init(XX, group, type, r,
                                  dropout, eps2, maxIter, output)
  if(!fast_version){
    # hH <- gfm1$hH; hB <- gfm1$hB; hmu <- gfm1$hmu
    gfm1 <-gfm_eval_intercept_osfinal(XX, gfm1$hH, gfm1$hB,gfm1$hmu, group, type)
  }
  hHm <- cbind(1, gfm1$hH)
  hBm <- cbind(gfm1$hmu, gfm1$hB)
  
  ICr <- sum(ICriteria(XX, hBm, hHm, r, group, type, criteria=cri))
  return(ICr)
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



signrevise <- function(A1, A2){
  nzid1 <- which(rowSums(A1^2)> 1e-5)[1]
  q <- ncol(A1)
  A <- sapply(1:q, function(k){
    if(sign(A1[nzid1,k]) != sign(A2[nzid1,k]))
      return(-A1[,k])
    return(A1[,k])
  })
  return(A)
}

Factorm <- function (X, q = NULL)
{
  if((!is.null(q)) && (q<1) ) stop("q must be NULL or other positive integer!")
  if(!is.matrix(X)) stop("X must be a matrix.")
  X <- scale(X, scale=FALSE)
  n <- nrow(X)
  p <- ncol(X)
  if (p > n) {
    svdX <- eigen(X %*% t(X))
    evalues <- svdX$values
    eigrt <- evalues[1:(21 - 1)]/evalues[2:21]
    if (is.null(q)) {
      q <- which.max(eigrt)
    }
    hatF <- as.matrix(svdX$vector[, 1:q] * sqrt(n))
    B2 <- n^(-1) * t(X) %*% hatF
    sB <- sign(B2[1, ])
    hB <- B2 * matrix(sB, nrow = p, ncol = q, byrow = TRUE)
    hH <- sapply(1:q, function(k) hatF[, k] * sign(B2[1,
    ])[k])
  }
  else {
    svdX <- eigen(t(X) %*% X)
    evalues <- svdX$values
    eigrt <- evalues[1:(21 - 1)]/evalues[2:21]
    if (is.null(q)) {
      q <- which.max(eigrt)
    }
    hB1 <- as.matrix(svdX$vector[, 1:q])
    hH1 <- n^(-1) * X %*% hB1
    svdH <- svd(hH1)
    hH2 <- signrevise(svdH$u * sqrt(n), hH1)
    if (q == 1) {
      hB1 <- hB1 %*% svdH$d[1:q] * sqrt(n)
    }
    else {
      hB1 <- hB1 %*% diag(svdH$d[1:q]) * sqrt(n)
    }
    sB <- sign(hB1[1, ])
    hB <- hB1 * matrix(sB, nrow = p, ncol = q, byrow = TRUE)
    hH <- sapply(1:q, function(j) hH2[, j] * sB[j])
  }
  sigma2vec <- colMeans((X - hH %*% t(hB))^2)
  res <- list()
  res$hH <- hH
  res$hB <- hB
  res$q <- q
  res$sigma2vec <- sigma2vec
  res$propvar <- sum(evalues[1:q])/sum(evalues)
  res$egvalues <- evalues
  attr(res, "class") <- "fac"
  return(res)
}

ortheB <- function(Z){
  Z=replace_with_column_means(Z)
  B <- qr(Z)
  eigvs <- sqrt(sort(eigen(t(Z)%*%Z)$values, decreasing = TRUE))
  B1 <- qr.Q(B) %*% Diag(eigvs)
  B0 <- B1 %*% Diag(sign(B1[1,]))
  return(B0)
}


ortheH <- function(H){
  H=replace_with_column_means(H)
  H1 <- qr.Q(qr(H)) * sqrt(nrow(H))
  hH <- H1 %*% Diag(sign(H[1,]) * sign(H1[1,]))
  return(hH)
}

family2func <- function(type){
  switch(type,
         gaussian= gaussian(link = "identity"),
         binomial = binomial(link = "logit"),
         poisson = poisson(link = "log"))
}


localupdateB2 <- function(X, g1, hH, type1){
  # g1 <- gcell[[1]]; type1 <- type[1]
  n <- nrow(X);q <- ncol(hH)
  p1 <- length(g1); B1 <- matrix(0, q+1, p1)
  jg <- 1:p1
  fun1 <- family2func(type1)
  #print(1)
  hH=replace_with_column_means(hH)
  #print(hH)
  hH[is.na(hH)] <- mean(hH, na.rm = TRUE)
  hH[is.infinite(hH)] <- NA  
  hH[is.na(hH)] <- mean(hH, na.rm = TRUE) 
  hH=t(apply(hH, 1, function(x) proj(x, 50)))
  
  for(j in jg){
    #hH=replace_with_column_means(hH)
    #print(j)
    B1[,j] <- glm.fit(x=cbind(1,hH), y=X[, g1[j]], family = fun1, intercept = FALSE)$coefficients
  }
  
  
  return(B1)
}

paraglmj1fit <- function(i,j, hBm, gcell, w, XX, funj){
  glm.fit(hBm[gcell[[j]], -1], XX[i, gcell[[j]]], weights = w,intercept = FALSE,
          family=funj,offset = hBm[gcell[[j]], 1])$coefficients
  
}

localupdateH2 <- function(X, gcell, hBm, type, dropout){
  n <- nrow(X); q <- ncol(hBm)-1; ng <- length(type)
  if(dropout !=0 && length(setdiff(dropout, 1:ng))>0){
    stop('dropout setting is wrong!')
  }
  idres <- setdiff(1:ng, dropout)
  Harray <- array(0, dim=c(n, q, length(idres)))
  w <- rep(1,n)
  for(jj in 1: length(idres)){
    j <- idres[jj]
    funj <- family2func(type[j])
    H2 <- matrix(rnorm(q*n), q, n)
    if(type[j] == 'gaussian'){
      w <- 1/ apply(X[,gcell[[j]]], 2, var)
      hBm=replace_with_column_means(hBm)
      for(i in 1:n){ #i <- 1
        H2[,i] <- glm.fit(hBm[gcell[[j]], -1], X[i, gcell[[j]]], weights = w,intercept = FALSE,
                          family=funj,offset = hBm[gcell[[j]], 1])$coefficients
      }
      
      
    }else{
      
      for(i in 1:n){
        try(
          H2[,i] <- glm.fit(hBm[gcell[[j]], -1], X[i, gcell[[j]]],intercept = FALSE,
                            family=funj,offset = hBm[gcell[[j]], 1])$coefficients
          , silent=TRUE
        )
      }
      
      
    }
    Harray[,,jj] <- t(H2)
  }
  hH <- apply(Harray, c(1,2), mean)
  hH=replace_with_column_means(hH)
  return(hH)
}



objfunc <- function(Hm, Bm, X, omega, gcell, type){
  n <- nrow(X); p <- ncol(X)
  eps1 <- 1e-20
  BHm <- Hm %*% t(Bm)
  ng <- length(type)
  Q <- matrix(0, n, p)
  for(j in 1:ng){
    if(type[j]== 'gaussian'){
      Q[,gcell[[j]]] <- (X[,gcell[[j]]] - BHm[,gcell[[j]]])^2
    }else if(type[j] == 'poisson'){
      me <- exp(BHm[,gcell[[j]]])
      Q[,gcell[[j]]] <- -log(dpois(X[, gcell[[j]]], me)+eps1)
    }else if(type[j] == 'binomial'){
      me3 <- 1 / (1 + exp(-BHm[,gcell[[j]]]))
      Q[,gcell[[j]]] <- -X[, gcell[[j]]] * log(me3+eps1) +
        (1-X[, gcell[[j]]]) * log(1-me3 + eps1)
    }
  }
  obj <- 1/n*omega*sum(Q)
  return(obj)
}


gfm_eval_intercept_init <- function(X, group, type, q,
                                    dropout, eps2, maxIter=10,
                                    output=0){
  ind_set <- unique(group)
  ng <- length(ind_set)
  if(length(setdiff(1:ng, ind_set))>0){
    stop("ID number of types must match type!")
  }
  if(ng != length(type)){
    stop("The number of groups must match with length of type!")
  }
  
  gcell <- list()
  for(j in 1:ng){
    gcell[[j]] <- which(group==j)
  }
  n <- nrow(X); p <- ncol(X)
  if(length(group) != p){
    stop("The length of group must match with column of X!")
  }
  omega <- 1/p
  #initialize
  hH <- Factorm(scale(X,scale=FALSE), q)$hH;hB <- 0
  eps1 <- 1e-4
  dBm <- Inf; dH <- Inf; dc =Inf; dOmega <- max(dBm, dH)
  tmpBm <- matrix(0, p, q+1); tmpH <- hH; tmpc <- 1e7
  k <- 1; history <- list()
  tic <- proc.time()
  while(k <= maxIter && dOmega > eps1 && dc >eps2){
    hhB <- NULL
    for(j in 1:ng){
      B1 <- localupdateB2(X, gcell[[j]], hH, type[j])
      hhB <- cbind(hhB, B1)
    }
    
    hmu <- hhB[1,]
    if(q == 1){
      hB <- matrix(hhB[-1,], ncol=1)
    }else{
      hB <- t(hhB[-1,])
    }
    hB=t(replace_with_column_means(t(hB)))
    #hmu=replace_with_column_means(hmu)
    # ensure indentifiability.
    hB <- ortheB(hB)
    hBm <- cbind(hmu, hB)
    dB <- norm(hBm-tmpBm, "F") / norm(hBm, 'F')
    tmpBm <- hBm
    if(output){
      message('---------- B updation is finished!---------\n')
    }
    H4 <- localupdateH2(X, gcell, hBm, type, dropout)
    H4=replace_with_column_means(H4)
    H4=t(replace_with_column_means(t(H4)))
    hH <- ortheH(H4)# %*% diag(sign(H4[1,]))
    dH <- norm(hH- tmpH, 'F')/norm(hH, 'F')
    tmpH <- hH
    if(output){
      message('---------- H updation is finished!---------\n')
    }
    hHm <- cbind(1, hH)
    hHm=replace_with_column_means(hHm)
    dOmega <- max(dB, dH)
    c <- objfunc(hHm, hBm, X, omega, gcell, type)
    dc <- abs(c - tmpc)/abs(tmpc)
    tmpc <- c
    if(output){
      message('Iter=', k, ', dB=',round(dB,4), ', dH=', round(dH,4),
              ',dc=', round(dc,4), ', c=', round(c,4), '\n')
    }
    history$dB[k] <- dB; history$dH[k] <- dH; history$dc[k] <- dc;
    history$c[k] <- c
    k <- k + 1
  }
  stoc <- proc.time() - tic
  history$realIter <- k -1
  history$maxIter <- maxIter
  history$elapsedTime <- stoc
  return(list(hH=hH, hB=hB, hmu=hmu, history=history))
}

gfm_eval_intercept_init_oneadd <- function(X, group, type, q,
                                           dropout, eps2, maxIter=10,
                                           output=0){
  ind_set <- unique(group)
  ng <- length(ind_set)
  if(length(setdiff(1:ng, ind_set))>0){
    stop("ID number of types must match type!")
  }
  if(ng != length(type)){
    stop("The number of groups must match with length of type!")
  }
  
  gcell <- list()
  for(j in 1:ng){
    gcell[[j]] <- which(group==j)
  }
  n <- nrow(X); p <- ncol(X)
  if(length(group) != p){
    stop("The length of group must match with column of X!")
  }
  omega <- 1/p
  #initialize
  hH <- Factorm(scale(X,scale=FALSE), q)$hH;hB <- 0
  eps1 <- 1e-4
  dBm <- Inf; dH <- Inf; dc =Inf; dOmega <- max(dBm, dH)
  tmpBm <- matrix(0, p, q+1); tmpH <- hH; tmpc <- 1e7
  k <- 1; history <- list()
  tic <- proc.time()
  while(k <= maxIter && dOmega > eps1 && dc >eps2){
    hhB <- NULL
    for(j in 1:ng){
      B1 <- localupdateB2(X, gcell[[j]], hH, type[j])
      hhB <- cbind(hhB, B1)
    }
    hhB=t(replace_with_column_means(t(hhB)))
    hmu <- hhB[1,]
    if(q == 1){
      hB <- matrix(hhB[-1,], ncol=1)
    }else{
      hB <- t(hhB[-1,])
    }
    # ensure indentifiability.
    #hB <- ortheB(hB)
    hBm <- cbind(hmu, hB)
    dB <- norm(hBm-tmpBm, "F") / norm(hBm, 'F')
    tmpBm <- hBm
    if(output){
      message('---------- B updation is finished!---------\n')
    }
    H4 <- localupdateH2(X, gcell, hBm, type, dropout)
    H4=replace_with_column_means(H4)
    #hH <- ortheH(H4)# %*% diag(sign(H4[1,]))
    hH <-H4
    dH <- norm(hH- tmpH, 'F')/norm(hH, 'F')
    tmpH <- hH
    if(output){
      message('---------- H updation is finished!---------\n')
    }
    hHm <- cbind(1, hH)
    dOmega <- max(dB, dH)
    c <- objfunc(hHm, hBm, X, omega, gcell, type)
    dc <- abs(c - tmpc)/abs(tmpc)
    tmpc <- c
    if(output){
      message('Iter=', k, ', dB=',round(dB,4), ', dH=', round(dH,4),
              ',dc=', round(dc,4), ', c=', round(c,4), '\n')
    }
    history$dB[k] <- dB; history$dH[k] <- dH; history$dc[k] <- dc;
    history$c[k] <- c
    k <- k + 1
  }
  stoc <- proc.time() - tic
  history$realIter <- k -1
  history$maxIter <- maxIter
  history$elapsedTime <- stoc
  return(list(hH=hH, hB=hB, hmu=hmu, history=history))
}



trace_statistic_fun <- function(H, H0){
  
  tr_fun <- function(x) sum(diag(x))
  mat1 <- t(H0) %*% H %*% ginv(t(H) %*% H) %*% t(H) %*% H0
  
  tr_fun(mat1) / tr_fun(t(H0) %*% H0)
  
}


measurefun <- function(hH, H, type=c('trace_statistic','ccor')){
  type <- match.arg(type)
  q <- min(ncol(H), ncol(hH))
  switch(type,
         ccor=cancor(hH, H)$cor[q],
         trace_statistic= trace_statistic_fun(hH, H))
}


localonestepH2 <- function(X, Bm, Hm, gcell, type){
  # Bm <- hBm_init ; Hm <- hHm_init
  n <- nrow(X); p <- ncol(X)
  if(ncol(Bm)-1 ==1){
    B <- matrix(Bm[,-1], ncol=1)
  }else{
    B <- Bm[,-1]
  }
  
  q <- ncol(Hm) - 1
  ng <- length(type)
  mucell <- list() # mean matrix
  for(j in 1:ng){
    mucell[[j]] <- switch(type[j],
                          gaussian = Hm %*% t(Bm[gcell[[j]],]),
                          poisson = exp(Hm %*% t(Bm[gcell[[j]],])),
                          binomial = 1/(1+exp(-Hm %*% t(Bm[gcell[[j]],])))
    )
  }
  # Score matrix
  df2 <- matrix(0,n, q)
  for(j in 1:ng){
    df2 <- df2 + (X[,gcell[[j]]] - mucell[[j]]) %*% B[gcell[[j]],]
  }
  H2update <- matrix(0, n, q)
  # Hessian matrix or information matrix
  # d2f <- list()
  for(i in 1:n){
    Bng <- matrix(0, q, q)
    for(j in 1:ng){
      dBng <- switch (type[j],
                      gaussian = t(B[gcell[[j]],]) %*% B[gcell[[j]],],
                      poisson = t(B[gcell[[j]],]) %*% diag(mucell[[j]][i,]) %*% B[gcell[[j]],],
                      binomial = t(B[gcell[[j]],]) %*% diag(mucell[[j]][i,]*(1-mucell[[j]][i,])) %*% B[gcell[[j]],]
      )
      Bng <- Bng + dBng
    }
    H2update[i,] <- qr.solve(Bng + 1e-6*diag(q)) %*% df2[i,]
  }
  H2 <- Hm[,-1] + H2update
  return(H2)
}

gfm_eval_intercept_osfinal <- function(X, hH_init, hB_init,hmu_init, group, type){
  # hH_init <- hH; hB_init <- hB; hmu_init <- hmu
  n <- nrow(hH_init); p <- nrow(hB_init)
  omega <- 1/p
  ind_set <- unique(group)
  ng <- length(ind_set)
  if(length(setdiff(1:ng, ind_set))>0){
    stop("ID number of types must match type!")
  }
  if(ng != length(type)){
    stop("The number of groups must match with length of type!")
  }
  if(ng == 1){
    hH_final <- hH_init;
    hB_final <- hB_init;
    hmu_final <- hmu_init;
    gcell <- list();gcell[[1]] <- 1:p
    objvalue <- objfunc(cbind(1, hH_final), cbind(hmu_final, hB_final),
                        X, omega, gcell, type)
  }else{
    gcell <- list()
    for(j in 1:ng){
      gcell[[j]] <- which(group ==j)
    }
    q <- ncol(hH_init)
    # update H
    hBm_init <- cbind(hmu_init, hB_init)
    hHm_init <- cbind(1, hH_init)
    hHm_init=replace_with_column_means(hHm_init)
    hHm_init=t(replace_with_column_means(t(hHm_init)))
    hBm_init=replace_with_column_means(hBm_init)
    hBm_init=t(replace_with_column_means(t(hBm_init)))
    H5 <- localonestepH2(X, hBm_init, hHm_init, gcell, type)
    H5=replace_with_column_means(H5)
    H5=t(replace_with_column_means(t(H5)))
    # update B
    hhB <- NULL
    for(j in 1:ng){
      B1 <- localupdateB2(X, gcell[[j]], H5, type[j])
      hhB <- cbind(hhB, B1)
    }
    hhB=t(replace_with_column_means(t(hhB)))
    hhB=replace_with_column_means(hhB)
    hmu <- hhB[1,]
    if(q == 1){
      hB <- matrix(hhB[-1,], ncol=1)
    }else{
      hB <- t(hhB[-1,])
    }
    
    # hBm_init <- cbind(hmu, hB)
    # hHm_init <- cbind(1, H5)
    # H5 <- localonestepH2(X, hBm_init, hHm_init, gcell, type)
    hH_final <- H5 # ortheH(H5)
    hB_final <- hB #ortheB(hB)
    hmu_final <- hmu
    objvalue <- objfunc(cbind(1, hH_final), cbind(hmu_final, hB_final),
                        X, omega, gcell, type)
  }
  return(list(hH=hH_final, hB=hB_final, hmu=hmu_final, obj=objvalue))
}


gfm_eval_intercept_osfinal_oneadd <- function(X, hH_init, hB_init,hmu_init, group, type){
  # hH_init <- hH; hB_init <- hB; hmu_init <- hmu
  n <- nrow(hH_init); p <- nrow(hB_init)
  omega <- 1/p
  ind_set <- unique(group)
  ng <- length(ind_set)
  if(length(setdiff(1:ng, ind_set))>0){
    stop("ID number of types must match type!")
  }
  if(ng != length(type)){
    stop("The number of groups must match with length of type!")
  }
  if(ng == 1){
    hH_final <- hH_init;
    hB_final <- hB_init;
    hmu_final <- hmu_init;
    gcell <- list();gcell[[1]] <- 1:p
    objvalue <- objfunc(cbind(1, hH_final), cbind(hmu_final, hB_final),
                        X, omega, gcell, type)
  }else{
    gcell <- list()
    for(j in 1:ng){
      gcell[[j]] <- which(group ==j)
    }
    q <- ncol(hH_init)
    # update H
    hBm_init <- cbind(hmu_init, hB_init)
    hHm_init <- cbind(1, hH_init)
    hBm_init=replace_with_column_means(hBm_init)
    hHm_init=replace_with_column_means(hHm_init)
    H5 <- localonestepH2(X, hBm_init, hHm_init, gcell, type)
    H5=replace_with_column_means(H5)
    # update B
    hhB <- NULL
    for(j in 1:ng){
      B1 <- localupdateB2(X, gcell[[j]], H5, type[j])
      hhB <- cbind(hhB, B1)
    }
    hmu <- hhB[1,]
    if(q == 1){
      hB <- matrix(hhB[-1,], ncol=1)
    }else{
      hB <- t(hhB[-1,])
    }
    
    # hBm_init <- cbind(hmu, hB)
    # hHm_init <- cbind(1, H5)
    # H5 <- localonestepH2(X, hBm_init, hHm_init, gcell, type)
    hB=replace_with_column_means(hB)
    hH_final <- H5 # ortheH(H5)
    hB_final <- hB #ortheB(hB)
    hmu_final <- hmu
    objvalue <- objfunc(cbind(1, hH_final), cbind(hmu_final, hB_final),
                        X, omega, gcell, type)
  }
  return(list(hH=hH_final, hB=hB_final, hmu=hmu_final, obj=objvalue))
}



chooseFacNumber_ratio <- function(XList, types, q_set = 2: 10, offset=FALSE, dc_eps=1e-4, maxIter=30,
                                  verbose = TRUE, parallelList=NULL){
  q_max <- max(q_set)
  gfm2 <- GFM::gfm(XList, types, q=q_max, algorithm = 'AM', offset=offset, dc_eps=dc_eps, maxIter=maxIter,verbose = verbose)
  svalues <- (svd(gfm2$hB))$d
  svalues_use <- svalues[svalues>1e-2]
  q_max <- length(svalues_use)
  q <- which.max(svalues[-q_max] / svalues[-1])
  message('The singular value ratio (SVR) method estimates the factor number q  as ', q, '. \n')
  return(q)
}



vem.fit <- function(XList, types, q, offset=FALSE, epsELBO=1e-5, maxIter=30, verbose=TRUE){
  # epsELBO=1e-5; maxIter=30; verbose=TRUE
  #typeID <- c(1, 2, 3);
  n <- nrow(XList[[1]]);
  pvec <- sapply(XList, ncol); p <- sum(pvec)
  A <- matrix(0, n, p);
  type_map <- 1:3;
  names(type_map) <- c("gaussian", "poisson", "binomial")
  typeID <- unname(type_map[types])
  
  if(offset && any(typeID==2)){
    
    AList <- XList
    A <- NULL
    for(i in seq_along(typeID)){
      dim1 <- dim(AList[[i]])
      AList[[i]] <- matrix(0, dim1[1], dim1[2])
      if(typeID[i] == 2){
        tmpvec <- rowSums(AList[[i]])
        AList[[i]] <- matrix(log(tmpvec+(tmpvec==0)), dim1[1], dim1[2])
      }
      A <- cbind(A, AList[[i]]) # re-obtain A.
    }
    
    rm(AList)
  }
  A <- as(A, "sparseMatrix")
  Mu_y_int <- NULL
  for(i in seq_along(typeID)){
    if(typeID[i]!=2){
      Mu_y_int <- cbind(Mu_y_int,XList[[i]])
    }else if(typeID[i] == 2){
      Mu_y_int <- cbind(Mu_y_int, log(1+ XList[[i]]))
    }
    
  }
  S_y_int = matrix(1, n, p);
  Fac <- get_initials(Mu_y_int, q= q)
  B_int <- Fac$hB
  Mu_h_int <- Fac$hH
  # B_int = matrix(rnorm(p*q), p, q)*0.1;
  # Mu_h_int <-  matrix(rnorm(n*q), n, q)#matrix(0, n, q)
  mu_int <- colMeans(Mu_y_int)
  #invLambda_int = rep(1, p);
  invLambda_int <- 1/ Fac$sigma2vec
  reslist <- VB_GFMcpp(XList, typeID, A,  Mu_y_int, S_y_int, invLambda_int,
                       B_int, mu_int, Mu_h_int, epsELBO, maxIter,
                       verbose)
  
  
  
  return(reslist)
}


VB_GFMcpp <- function(XList, typeID, A, Mu_y_int, S_y_int, invLambda_int, B_int, mu_int, H_int, epsELBO, maxIter, verbose) {
  .Call(`_GFM_VB_GFMcpp`, XList, typeID, A, Mu_y_int, S_y_int, invLambda_int, B_int, mu_int, H_int, epsELBO, maxIter, verbose)
}



get_initials <- function(X, q){
  
  n <- nrow(X); p <- ncol(X)
  mu <- colMeans(X)
  X <- X - matrix(mu, nrow=n, ncol=p, byrow=TRUE)
  svdX  <- irlba(A =X, nv = q)
  PCs <- sqrt(n) * svdX$u
  loadings <- svdX$v %*% Diag(svdX$d[1:q]) / sqrt(n)
  dX <- PCs %*% t(loadings) - X
  Lam_vec <- colSums(dX^2)/n
  return(list(hH = PCs, hB = loadings, hmu=mu,sigma2vec = Lam_vec))
  
}
library(Rcpp)
library(RcppArmadillo)
Rcpp::sourceCpp('~/vsdt/code/realdata/bgfm.cpp')


# obtain the value of IC4,IC5,ICDT,VSDT
fac_ic <- function(x, qmax = 20, reptimes = 10, maxit=30) {
  fis_mat = c()
  cov_mat <- matrix(NA, ncol = qmax, nrow = reptimes)
  loglik_mat = matrix(NA, ncol = qmax, nrow = reptimes)
  n = dim(x)[1]
  p = dim(x)[2]
  
  for (b in 1:reptimes) {
    set.seed(b + 1000)
    
    
    XA = matrix(nrow = n, ncol = p)
    XB = matrix(nrow = n, ncol = p)
    
    for (i in 1:n) {
      for (j in 1:p) {
        XB[i, j] = x[i, j] - rbinom(1, x[i, j], 0.5)
      }
    }
    
    XA = x - XB
    lnlikB = vector(mode = "numeric", length = qmax)
    temp_cov_mat = vector(mode = "numeric", length = qmax)
    
    tryCatch({
      for (r in 1:qmax) {
        
        
        re0326 <- GFM::gfm(list(XA), types = c("poisson"), q = r, algorithm = "AM",maxIter = maxit, verbose = F)
        para <- exp(re0326$hH %*% t(re0326$hB) + rep(1, n) %*% t(re0326$hmu))
        lnlikB[r] = sum(-log(dpois(XB, para) + 1e-20)) / (n * p)
        temp_cov_mat[r] = lnlikB[r] - re0326$obj
        #cat("r = ", r,"\n")
      }
    }, error = function(e) {
      
      for (r in 1:qmax) {
        re0326 <<- gfm(list(XA), types = c("poisson"), q = r,maxIter = maxit, verbose = F)
        para <<- exp(re0326$hH %*% t(re0326$hB) + rep(1, n) %*% t(re0326$hmu))
        lnlikB[r] <<- sum(-log(dpois(XB, para) + 1e-20)) / (n * p)
        temp_cov_mat[r] <<- lnlikB[r] - re0326$obj
        #cat("r = ", r,"\n")
      }
    })
    
    cov_mat[b,] = temp_cov_mat
    fis_mat[b] = which.min(lnlikB)
    loglik_mat[b,] = lnlikB
    cat("VSDT result = ", fis_mat[b], "\n")
    # cat("lnlik B = ",lnlikB,"\n")
    # cat("b = ", b,"\n")
  }
  
  cat("VSDT \n")
  print(fis_mat)
  log_VSDT = loglik_mat[1,]
  cat("VSDT value = \n")
  print(loglik_mat)
  
  cov_vec = colMeans(cov_mat,na.rm = T)
  log_IC4 = c()
  log_IC5 = c()
  log_DT = c()
  logli = c()
  
  n = dim(x)[1]
  p = dim(x)[2]
  
  tryCatch({
    for (r in 1:qmax) {
      
      
      re0326_all = GFM::gfm(list(x), types = c("poisson"), q = r, algorithm = "AM", maxIter = maxit, verbose = F)
      logli[r] = re0326_all$obj
      log_IC5[r] = re0326_all$obj + r * ((n + p) / (n * p)) * log(n * p / (n + p))
      log_IC4[r] = re0326_all$obj + r / min(sqrt(n), sqrt(p))^2 * log(min(sqrt(n), sqrt(p))^2)
      log_DT[r] = re0326_all$obj + cov_vec[r] * log(n * p / (n + p))
      cat("r = ", r,"\n")
    }
  }, error = function(e) {
    #cat("Error encountered in the second part, switching to the backup method.\n")
    
    for (r in 1:qmax) {
      re0326_all <<- gfm(list(x), types = c("poisson"), q = r, maxIter = maxit, verbose = F)
      logli[r] <<- re0326_all$obj
      log_IC5[r] <<- re0326_all$obj + r * ((n + p) / (n * p)) * log(n * p / (n + p))
      log_IC4[r] <<- re0326_all$obj + r / min(sqrt(n), sqrt(p))^2 * log(min(sqrt(n), sqrt(p))^2)
      log_DT[r] <<- re0326_all$obj + cov_vec[r] * log(n * p / (n + p))
      #cat("r = ", r,"\n")
    }
  })
  
  return(list(
    loglike_val = logli,
    IC5_val = log_IC5,
    IC4_val = log_IC4,
    ICDT_val = log_DT,
    VSDT_val_mat = loglik_mat
  ))
}





iden5q0309=function(H,B,c,mu,y){
  q = ncol(H); n = nrow(H)
  mu525=mu+(B)%*%colMeans(H)
  svdHB525 <- irlba((H- matrix(colMeans(H), n, 5, byrow = TRUE)) %*% t(B), nv= 5)
  signB1525 <-c(apply(svdHB525$v, 2, function(x) sign(x[which.max(abs(x))])))
  H525 <- sqrt(n) * svdHB525$u %*% diag(c(signB1525))
  B525 <- (svdHB525$v %*% diag(c(svdHB525$d*signB1525)) / sqrt(n))
  lmw1=lm(H525[,1]~y-1)
  lmw2=lm(H525[,2]~y-1)
  lmw3=lm(H525[,3]~y-1)
  lmw4=lm(H525[,4]~y-1)
  lmw5=lm(H525[,5]~y-1)
  CC525=c+B525[,1]*(lmw1$coefficients)+B525[,2]*(lmw2$coefficients)+B525[,3]*(lmw3$coefficients)+B525[,4]*(lmw4$coefficients)+B525[,5]*(lmw5$coefficients)
  HH525=cbind(lmw1$residuals,lmw2$residuals,lmw3$residuals,lmw4$residuals,lmw5$residuals)
  BB525=B525
  mu525=mu525
  return(list(H=HH525, B=BB525,c=CC525, mu=mu525))
}


# present results
calculate_proportions <- function(x) {
  valid_x <- x[x != 0 & !is.na(x)]  
  total_count <- length(valid_x)  
  
  rounded_equal_5 <- round(sum(valid_x == 5) / total_count * 100)
  rounded_between_1_and_5 <- round(sum(valid_x >= 1 & valid_x < 5) / total_count * 100)
  rounded_greater_than_5 <- 100 - (rounded_equal_5 + rounded_between_1_and_5)
  
  c("Correct estimation" = rounded_equal_5, "Over estimation" = rounded_greater_than_5, "Under estimation" = rounded_between_1_and_5)
}

