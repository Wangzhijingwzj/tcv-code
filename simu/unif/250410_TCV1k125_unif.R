

library(foreach)
library(doParallel)
library(irlba)
library(doSNOW)
library(openxlsx)
library(countsplit)
library(bcv)
library(tidyr)
library(invgamma)
library(mvtnorm)
library(pracma)
library(pcaPP)
library(Rfast)
library(dplyr)
# Create cluster and register parallel backend
cl <- makeCluster(40)
registerDoParallel(cl)
#RhpcBLASctl::blas_set_num_threads(1)


# Function to extract factor numbers
extract_factor_numbers <- function(results_list) {
  # Initialize an empty matrix to store factor number estimates
  facnum_mat <- matrix(NA, nrow = length(results_list), 
                       ncol = 33)  # 32 columns for different criteria + 1 for ratio
  
  # Column names for clarity
  colnames(facnum_mat) <- c(
    "facnum_ratio", "facnum_bcv", "facnum_bcv2","facnum_sr",
    "facnum_IC1", "facnum_IC2", "facnum_zero", 
    "facnum_AIC1", "facnum_AIC2", "facnum_AIC3", 
    "facnum_BIC1", "facnum_BIC2", "facnum_BIC3",
    "facnum_IIC1", "facnum_IIC2", 
    "facnum_IAIC1", "facnum_IAIC2", "facnum_IAIC3",
    "facnum_IBIC1", "facnum_IBIC2", "facnum_IBIC3",
    "facnum_TCV2K", "facnum_TPCV2K", "facnum_TICV2K", "facnum_TIPCV2K",
    "facnum_TCV5K", "facnum_TPCV5K", "facnum_TICV5K", "facnum_TIPCV5K",
    "facnum_TCV10K", "facnum_TPCV10K", "facnum_TICV10K", "facnum_TIPCV10K"
  )
  
  # Iterate through each iteration in the results list
  for (s in seq_along(results_list)) {
    # Extract the result matrix for this iteration
    result_matrix <- results_list[[s]]$result_matrix
    
    # 1. Extract factor numbers based on various information criteria
    info_criteria <- c(
      "log_IC1", "log_IC2", "log_zero", 
      "log_AIC1", "log_AIC2", "log_AIC3", 
      "log_BIC1", "log_BIC2", "log_BIC3",
      "log_IIC1", "log_IIC2", 
      "log_IAIC1", "log_IAIC2", "log_IAIC3",
      "log_IBIC1", "log_IBIC2", "log_IBIC3"
    )
    
    # Calculate factor numbers for each information criterion
    for (i in seq_along(info_criteria)) {
      facnum_mat[s, i+4] <- which.min(result_matrix[, info_criteria[i]])
    }
    
    # 2. Extract other factor number estimates
    facnum_mat[s, 1] <- results_list[[s]]$facnum_ratio
    facnum_mat[s, 2] <- results_list[[s]]$facnum_bcv
    facnum_mat[s, 3] <- results_list[[s]]$facnum_bcv2
    facnum_mat[s, 4] <- results_list[[s]]$facnum_sr
    # 3. Extract factor numbers based on cross-validation metrics
    cv_metrics <- c(
      "log_TCV2K", "log_TPCV2K", "log_TICV2K", "log_TIPCV2K",
      "log_TCV5K", "log_TPCV5K", "log_TICV5K", "log_TIPCV5K",
      "log_TCV10K", "log_TPCV10K", "log_TICV10K", "log_TIPCV10K"
    )
    
    # Calculate factor numbers for each cross-validation metric
    for (i in seq_along(cv_metrics)) {
      facnum_mat[s, i+21] <- which.min(results_list[[s]][[cv_metrics[i]]])
    }
  }
  
  return(facnum_mat)
}

# Utility Functions
getModeAverage <- function(v) {
  # Calculate mode or average of mode for a vector
  tab <- table(v)
  modes <- names(tab[tab == max(tab)])
  
  if (length(modes) > 1) {
    # Multiple modes: calculate average
    modes <- as.numeric(modes)
    modeAverage <- mean(modes)
  } else {
    # Single mode: convert to numeric if original vector is numeric
    if (is.numeric(v)) modes <- as.numeric(modes)
    modeAverage <- modes
  }
  
  return(modeAverage)
}

has_full_zero_row_or_col <- function(mat) {
  # Check if matrix has any full zero rows or columns
  any(rowSums(mat == 0) == ncol(mat)) || any(colSums(mat == 0) == nrow(mat))
}
cat("unif s=100, n1000 p100,200,500")


# Main Simulation Loop
S <- 100  # Number of iterations
for (n in c(1000)) {
  for (p in c(100,200,500)) {
    # Print current simulation parameters
    print("unif")
    print(paste("n:", n))
    print(paste("p:", p))
    
    # Parallel computation
    results_list <- foreach(
      s = 1:S, 
      .packages = c("irlba", "doSNOW", "countsplit", "foreach", "doParallel", "GFM","bcv")
    ) %dopar% {
      # Source additional functions
      source("~/multiDT/250113function_tcv.R")
      Rcpp::sourceCpp("~/glmcv/simu/confirm_jmle_omp_poisson_with_intercept_missing.cpp")
      
      # Create a list to store all results for this iteration
      all_results <- list()
      
      # Set seed for reproducibility
      set.seed(s + 1000)
      
      # Factor matrix generation
      q <- 5
      FF <- matrix(runif(n * q, min = -2, max = 2), nrow = n, ncol = q)
      all_results$q <- q
      
      
      # Loading matrix generation
      BB <- matrix(runif(p * q, min = -1, max = 1), nrow = p, ncol = q)
      
      
      # Additional parameter generation
      a <- runif(p, min = 0, max = 1)
      
      
      # Identifiability adjustment
      FF0 <- add_identifiability(FF, BB, a)$H
      BB0 <- add_identifiability(FF, BB, a)$B
      alpha <- add_identifiability(FF, BB, a)$mu
      
      
      # Lambda parameter calculation
      lambda <- exp(FF0 %*% t(BB0) + as.matrix(rep(1, n)) %*% t(as.matrix(alpha)))
      
      
      #perb n
      #lambda[((0.9*n+1):n),]=2*lambda[((0.9*n+1):n),]
      #perb p
      #lambda[,((0.9*p+1):p)]=2*lambda[,((0.9*p+1):p)]
      
      # Data generation
      x <- matrix(nrow = n, ncol = p)
      for (i in 1:n) {
        for (j in 1:p) {
          x[i, j] <- rpois(1, lambda = lambda[i, j])
        }
      }
      all_results$FF <- FF
      all_results$BB <- BB
      all_results$a <- a
      all_results$FF0 <- FF0
      all_results$BB0 <- BB0
      all_results$alpha <- alpha
      all_results$lambda <- lambda
      all_results$x <- x
      
      ###### ratio
      facnum_ratio <- tryCatch({
        suppressMessages(
          chooseFacNumber_ratio(list(x), types = c("poisson"), 
                                q_set = 1:8, select_method = "SVR", verbose = FALSE)
        )
      }, error = function(e){
        NA
      })
      all_results$facnum_ratio <- facnum_ratio
      
      
      #####cv1 
      facnum_bcv <- tryCatch({
        re_bcv = cv.svd.gabriel(sweep(log(x + 0.5), 1, rowMeans(log(x + 0.5)), "-"), krow = 2, kcol = 2, maxrank = 8)
        as.numeric(names(which.min(colMeans(re_bcv$msep)[-1])))
      }, error = function(e) {
        
        cat("Error occurred in cv1 (Gabriel method):\n")
        print(e)
        if(exists("s")) {
          cat("Current s value:", s, "\n")
        }
        
        NA
      })
      
      
      all_results$facnum_bcv <- facnum_bcv
      
      ###### cv2
      facnum_bcv2 <- tryCatch({
        re_bcv2 = cv.svd.wold(sweep(log(x + 0.5), 1, rowMeans(log(x + 0.5)), "-"), k = 5, maxrank = 8, tol = 1e-03, maxiter = 20)
        as.numeric(names(which.min(colMeans(re_bcv2$msep)[-1])))
      }, error = function(e) {
        
        cat("Error occurred in cv2 (Wold method):\n")
        print(e)
        if(exists("s")) {
          cat("Current s value:", s, "\n")
        }
        
        NA
      })
      
      all_results$facnum_bcv2 <- facnum_bcv2
      
      
      
      ######### sr 
      facnum_sr <- tryCatch({
        Rho <- cor(sweep(log(x + 0.5), 1, rowMeans(log(x + 0.5)), "-"), method = "spearman")
        eigvals_Rho <- eig(Rho)
        SR(eigvals_Rho, Kmax=8, p, n)
      }, error = function(e) {
        
        cat("Error occurred in sr (qiu 2024 jasa):\n")
        print(e)
        if(exists("s")) {
          cat("Current s value:", s, "\n")
        }
        
        NA
      })
      
      all_results$facnum_sr <- facnum_sr
      
      ###### Compute various cross-validation metrics
      for (k in c(2, 5, 10)) {
        tryCatch({
          re_TCV <- multiDT(x, K = k, rmax = 8)
          
          # Store results for each k
          all_results[[paste0("log_TCV", k, "K")]] <- re_TCV$TCV
          all_results[[paste0("log_TPCV", k, "K")]] <- 
            re_TCV$TCV + seq(1:8) * ((n+p)/(n*p)) * log((n*p)/(n+p))
          all_results[[paste0("log_TICV", k, "K")]] <- re_TCV$TICV
          all_results[[paste0("log_TIPCV", k, "K")]] <- 
            re_TCV$TICV + seq(1:8) * ((n+p)/(n*p)) * log((n*p)/(n+p))
        }, error = function(e) {
          for (metric in c("TCV", "TPCV", "TICV", "TIPCV")) {
            all_results[[paste0("log_", metric, k, "K")]] <- rep(NA, 8)
          }
        })
      }
      
      # Parallel computation of factor number estimates
      results <- foreach(r = 1:8) %dopar% {
        source("~/multiDT/250113function_tcv.R")
        Rcpp::sourceCpp("~/glmcv/simu/confirm_jmle_omp_poisson_with_intercept_missing.cpp")
        
        tryCatch({
          re0326_all <- tryCatch({
            GFM::gfm(list(x), types = c("poisson"), maxIter = 30, q = r, algorithm = "AM", verbose = FALSE)
          }, error = function(e) {
            confirm_CJMLE_poisson_cpp(x, matrix(TRUE, nrow = n, ncol = p), 
                                      cbind(rep(1, n), matrix(runif(n * r, 0, 1), n, r)), 
                                      matrix(runif(p * (r + 1), -1, 1), p, (r + 1)), 
                                      matrix(TRUE, p, (r + 1)), C = 100, tol = 0.01)
          })
          
          # Compute various information criteria
          list(
            log_IC1 = 2*re0326_all$obj + r * ((n + p) / (n * p)) * log(n * p / (n + p)),
            log_IC2 = 2*re0326_all$obj + r / min(sqrt(n), sqrt(p))^2 * log(min(sqrt(n), sqrt(p))^2),
            log_zero = 2*re0326_all$obj,
            log_AIC1 = 2*re0326_all$obj + r * 2 / n,
            log_AIC2 = 2*re0326_all$obj + r * 2 / p,
            log_AIC3 = 2*re0326_all$obj + r * 2 * (n + p - r) / (n * p),
            log_BIC1 = 2*re0326_all$obj + r * log(n) / n,
            log_BIC2 = 2*re0326_all$obj + r * log(p) / p,
            log_BIC3 = 2*re0326_all$obj + r * (n + p - r) * log(n * p) / (n * p),
            log_IIC1 = 2*log(re0326_all$obj) + r * ((n + p) / (n * p)) * log(n * p / (n + p)),
            log_IIC2 = 2*log(re0326_all$obj) + r / min(sqrt(n), sqrt(p))^2 * log(min(sqrt(n), sqrt(p))^2),
            log_IAIC1 = 2*log(re0326_all$obj) + r * 2 / n,
            log_IAIC2 = 2*log(re0326_all$obj) + r * 2 / p,
            log_IAIC3 = log(2*re0326_all$obj) + r * 2 * (n + p - r) / (n * p),
            log_IBIC1 = log(2*re0326_all$obj) + r * log(n) / n,
            log_IBIC2 = log(2*re0326_all$obj) + r * log(p) / p,
            log_IBIC3 = log(2*re0326_all$obj) + r * (n + p - r) * log(n * p) / (n * p),
            re0326_all = re0326_all
          )
        }, error = function(e) {
          # Return NA values if computation fails
          lapply(result_names, function(x) NA)
        })
      }
      
      # Extract results from the list
      result_matrix <- do.call(rbind, lapply(results, function(x) {
        if(is.list(x)) {
          x[names(x) != "re0326_all"]
        } else {
          x
        }
      }))
      
      # Store re0326_all results for each r
      all_results$re0326_all_results <- lapply(results, function(x) x$re0326_all)
      
      # Combine all results for this iteration
      all_results$result_matrix <- result_matrix
      
      return(all_results)
    }
    
    # Process factor numbers
    facnum_mat <- extract_factor_numbers(results_list)
    
    # Convert the list into a data frame
    res <- as.data.frame(facnum_mat)
    
    # Count occurrences of values 1 to 8 in each column of the data frame
    count_data <- sapply(res, function(x) table(factor(x, levels = 1:8)))
    
    # Rename columns to match the original data frame names
    colnames(count_data) <- names(res)
    
    # Combine the count data with the original data
    final_data <- rbind(count_data, res)
    
    # Prepare output directory
    setwd("~/multiDT/250409_simu/unif")
    
    # Create the raw results directory if it doesn't exist
    dir.create("raw_results", showWarnings = FALSE)
    
    # File names
    file_name_excel <- paste0("funif_n", n, "_p", p, "_q5_with_counts.xlsx")
    file_name_rds_counts <- paste0("funif_n", n, "_p", p, "_q5_with_counts.rds")
    file_name_rds_raw <- paste0("raw_results/funif_n", n, "_p", p, "_q5_raw_results.rds")
    
    # Create Excel workbook and save
    wb <- createWorkbook()
    addWorksheet(wb, "Sheet 1")
    writeData(wb, "Sheet 1", final_data)
    saveWorkbook(wb, file_name_excel, overwrite = TRUE)
    
    # Save count results as RDS
    saveRDS(final_data, file_name_rds_counts)
    print(final_data)
    # Save raw results as RDS
    saveRDS(results_list, file_name_rds_raw)
    
    cat("Files saved:\n")
    cat("Excel:", file_name_excel, "\n")
    cat("Counts RDS:", file_name_rds_counts, "\n")
    cat("Raw Results RDS:", file_name_rds_raw, "\n")
  }
}

# Clean up parallel cluster
stopCluster(cl)















