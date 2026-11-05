# # Load required libraries
# library(foreach)
# library(doParallel)
# library(irlba)
# library(doSNOW)
# library(openxlsx)
# library(countsplit)
# library(bcv)
# library(tidyr)
# library(invgamma)
# library(mvtnorm)
# library(pracma)
# library(pcaPP)
# library(Rfast)
# library(dplyr)
# library(ICtest)
# 
# # Create cluster and register parallel backend if needed
# # cl <- makeCluster(20)
# # registerDoParallel(cl)
# # RhpcBLASctl::blas_set_num_threads(1)
# 
# has_full_zero_row_or_col <- function(mat) {
#   # Check if matrix has any full zero rows or columns
#   any(rowSums(mat == 0) == ncol(mat)) || any(colSums(mat == 0) == nrow(mat))
# }
# 
# # Function to run the simulation
# run_simulation <- function(data_gen_method = "perbn", sample_sizes = c(100, 500), dimensions = c(100, 500), num_iterations = 10) {
#   # Create base output directory
#   base_output_dir <- "~/multiDT/250409_simu"
#   
#   cat("Running simulation for", data_gen_method, "with s =", num_iterations, 
#       ", n =", paste(sample_sizes, collapse=", "), 
#       ", p =", paste(dimensions, collapse=", "), "\n")
#   
#   # Main Simulation Loop
#   S <- num_iterations
#   for (n in sample_sizes) {
#     for (p in dimensions) {
#       # Print current simulation parameters
#       cat("\n=============================\n")
#       cat("Scenario:", data_gen_method, "\n")
#       cat("Sample size (n):", n, "\n")
#       cat("Dimension (p):", p, "\n")
#       cat("=============================\n")
#       
#       # Create storage for all results
#       all_results <- list()
#       facnum_act_all <- c()
#       facnum_ladel_all <- c()
#       
#       for (s in 1:S) {
#         cat("Iteration", s, "of", S, "\n")
#         
#         # Source function file
#         source("~/multiDT/250113function_tcv.R")
#         
#         # Set seed for reproducibility
#         set.seed(s + 1000)
#         
#         # Factor matrix generation based on method
#         if (data_gen_method == "q3noise" || data_gen_method == "q2noise") {
#           # For q3noise and q2noise methods, set q accordingly
#           q <- if(data_gen_method == "q3noise") 3 else 2
#           FF <- matrix(rnorm(n * q), nrow = n, ncol = q)
#         } else if (data_gen_method == "unif") {
#           # For unif method, use uniform distribution for FF
#           q <- 5
#           FF <- matrix(runif(n * q, min = -2, max = 2), nrow = n, ncol = q)
#         } else {
#           # Default: use normal distribution with q=5
#           q <- 5
#           FF <- matrix(rnorm(n * q), nrow = n, ncol = q)
#         }
#         
#         # Store q value in all_results
#         all_results$q <- q
#         
#         # Loading matrix generation (same for all methods)
#         BB <- matrix(runif(p * q, min = -1, max = 1), nrow = p, ncol = q)
#         
#         # Additional parameter generation (same for all methods)
#         a <- runif(p, min = 0, max = 1)
#         
#         # Identifiability adjustment
#         result <- add_identifiability(FF, BB, a)
#         FF0 <- result$H
#         BB0 <- result$B
#         alpha <- result$mu
#         
#         # Lambda parameter calculation based on method
#         if (data_gen_method == "q3noise") {
#           # Add noise to lambda for q3noise
#           lambda <- exp(FF0 %*% t(BB0) + as.matrix(rep(1, n)) %*% t(as.matrix(alpha)) + 
#                           matrix(rnorm(n * p, mean = 0, sd = 0.3), nrow = n, ncol = p))
#         } else if (data_gen_method == "q2noise") {
#           # Add noise to lambda for q2noise
#           lambda <- exp(FF0 %*% t(BB0) + as.matrix(rep(1, n)) %*% t(as.matrix(alpha)) + 
#                           matrix(rnorm(n * p, mean = 0, sd = 0.2), nrow = n, ncol = p))
#         } else if (data_gen_method == "norm") {
#           # Add noise to lambda for norm
#           lambda <- exp(FF0 %*% t(BB0) + as.matrix(rep(1, n)) %*% t(as.matrix(alpha)) + 
#                           matrix(rnorm(n * p, mean = 0, sd = 0.5), nrow = n, ncol = p))
#         } else {
#           # Default lambda calculation for other methods
#           lambda <- exp(FF0 %*% t(BB0) + as.matrix(rep(1, n)) %*% t(as.matrix(alpha)))
#         }
#         
#         # Store lambda in all_results
#         all_results$lambda <- lambda
#         
#         # Apply specific modifications based on method
#         if (data_gen_method == "perbn") {
#           # Perturb n (10% of the observations)
#           lambda[((0.9*n+1):n),] <- 2*lambda[((0.9*n+1):n),]
#         } else if (data_gen_method == "perbp") {
#           # Perturb p (10% of the variables)
#           lambda[,((0.9*p+1):p)] <- 2*lambda[,((0.9*p+1):p)]
#         }
#         
#         # Data generation
#         x <- matrix(nrow = n, ncol = p)
#         
#         if (data_gen_method == "over") {
#           # Mixed distribution (overdispersion) data generation
#           random_matrix <- runif(n*p)
#           
#           selection_threshold1 <- quantile(random_matrix, 0.1)
#           selection_threshold2 <- quantile(random_matrix, 0.2)
#           selection_threshold3 <- quantile(random_matrix, 0.3)
#           
#           for (i in 1:n) {
#             for (j in 1:p) {
#               if (random_matrix[(i-1)*p + j] <= selection_threshold1) {
#                 # For random values <= threshold1, use theta = 0.1
#                 theta <- 0.1
#                 size <- theta
#                 prob <- size / (size + lambda[i, j])
#                 x[i, j] <- rnbinom(1, size, prob)
#               } else if ((random_matrix[(i-1)*p + j] <= selection_threshold2) & 
#                          (random_matrix[(i-1)*p + j] > selection_threshold1)) {
#                 # For random values <= threshold2 but > threshold1, use theta = 1
#                 theta <- 1
#                 size <- theta
#                 prob <- size / (size + lambda[i, j])
#                 x[i, j] <- rnbinom(1, size, prob)
#               } else if ((random_matrix[(i-1)*p + j] <= selection_threshold3) & 
#                          (random_matrix[(i-1)*p + j] > selection_threshold2)) {
#                 # For random values <= threshold3 but > threshold2, use theta = 10
#                 theta <- 10
#                 size <- theta
#                 prob <- size / (size + lambda[i, j])
#                 x[i, j] <- rnbinom(1, size, prob)
#               } else {
#                 # For random values > threshold3, generate Poisson distribution result
#                 x[i, j] <- rpois(1, lambda[i, j])
#               }
#             }
#           }
#         } else {
#           # Standard Poisson data generation for other methods
#           for (i in 1:n) {
#             for (j in 1:p) {
#               x[i, j] <- rpois(1, lambda = lambda[i, j])
#             }
#           }
#         }
#         
#         # Apply zero-inflation for zip03 method
#         if (data_gen_method == "zip03") {
#           repeat {
#             # Create a copy of the original matrix to preserve non-zero structure
#             x_temp <- x
#             
#             # Calculate number of elements to set to zero
#             num_zeros <- ceiling(0.3 * length(x_temp))  # 30% of elements
#             
#             # Randomly select these element indices
#             indices <- sample(length(x_temp), num_zeros, replace = FALSE)
#             
#             # Set selected elements to 0
#             x_temp[indices] <- 0
#             
#             # Check if there are any full zero rows or columns
#             if (!has_full_zero_row_or_col(x_temp)) {
#               x <- x_temp  # Update x to new matrix
#               break        # Exit loop
#             }
#           }
#         }
#         
#         # Store generated data
#         all_results[[s]] <- list(
#           q = q,
#           FF = FF,
#           BB = BB,
#           a = a,
#           FF0 = FF0,
#           BB0 = BB0,
#           alpha = alpha,
#           lambda = lambda,
#           x = x
#         )
#         
#         # Preprocess data - centering log-transformed data
#         log_centered_x <- sweep(log(x + 0.5), 1, rowMeans(log(x + 0.5)), "-")
#         
#         # ACT method for factor number estimation
#         facnum_act <- tryCatch({
#           R <- cor(log_centered_x, method = "pearson")
#           eigvals_R <- eig(R)
#           ACT(eigvals_R, (1 + sqrt(p / n)), Kmax=8)
#         }, error = function(e) {
#           cat("Error occurred in ACT (Fan 2021 JASA):\n")
#           print(e)
#           cat("Current iteration:", s, "\n")
#           NA
#         })
#         
#         # LADEL method for factor number estimation
#         facnum_ladel <- tryCatch({
#           PCAladle(log_centered_x, ncomp = 8)$k
#         }, error = function(e) {
#           cat("Error occurred in LADEL (Luo 2016):\n")
#           print(e)
#           cat("Current iteration:", s, "\n")
#           NA
#         })
#         
#         # Store factor number estimates
#         facnum_act_all[s] <- facnum_act
#         facnum_ladel_all[s] <- facnum_ladel
#         
#         all_results[[s]]$facnum_act <- facnum_act
#         all_results[[s]]$facnum_ladel <- facnum_ladel
#       }
#       
#       # Create results matrix
#       facnum_mat <- cbind(
#         facnum_act = facnum_act_all,
#         facnum_ladel = facnum_ladel_all
#       )
#       
#       # Convert to data frame
#       res_df <- as.data.frame(facnum_mat)
#       
#       # Count occurrences of values 1 to 8 in each column
#       count_data <- sapply(res_df, function(x) {
#         tab <- table(factor(x, levels = 1:8))
#         if (length(tab) < 8) {
#           tab_full <- rep(0, 8)
#           names(tab_full) <- as.character(1:8)
#           tab_full[names(tab)] <- tab
#           return(tab_full)
#         } else {
#           return(tab)
#         }
#       })
#       
#       # Combine the count data with the original results
#       final_data <- rbind(count_data, res_df)
#       
#       # Prepare output directory and file names based on data generation method
#       if (data_gen_method == "perbn") {
#         output_dir <- "~/multiDT/250409_simu/perbn"
#         dir.create(output_dir, showWarnings = FALSE)
#         dir.create(file.path(output_dir, "raw_results"), showWarnings = FALSE)
#         
#         file_prefix <- paste0("act_ladel_fnorm_perbn_n", n, "_p", p, "_q5")
#         file_name_excel <- file.path(output_dir, paste0(file_prefix, "_with_counts.xlsx"))
#         file_name_rds_counts <- file.path(output_dir, paste0(file_prefix, "_with_counts.rds"))
#         file_name_rds_raw <- file.path(output_dir, "raw_results", paste0(file_prefix, "_raw_results.rds"))
#       } else if (data_gen_method == "perbp") {
#         output_dir <- "~/multiDT/250409_simu/perbp"
#         dir.create(output_dir, showWarnings = FALSE)
#         dir.create(file.path(output_dir, "raw_results"), showWarnings = FALSE)
#         
#         file_name_excel <- file.path(output_dir, paste0("act_ladel_fnorm_perbp_n", n, "_p", p, "_q5_with_counts.xlsx"))
#         file_name_rds_counts <- file.path(output_dir, paste0("act_ladel_fnorm_perbp_n", n, "_p", p, "_q5_with_counts.rds"))
#         file_name_rds_raw <- file.path(output_dir, paste0("raw_results/act_ladel_fnorm_perbp_n", n, "_p", p, "_q5_raw_results.rds"))
#       } else if (data_gen_method == "over") {
#         output_dir <- "~/multiDT/250409_simu/over"
#         dir.create(output_dir, showWarnings = FALSE)
#         dir.create(file.path(output_dir, "raw_results"), showWarnings = FALSE)
#         
#         file_name_excel <- file.path(output_dir, paste0("act_ladel_fnorm_over_n", n, "_p", p, "_q5_with_counts.xlsx"))
#         file_name_rds_counts <- file.path(output_dir, paste0("act_ladel_fnorm_over_n", n, "_p", p, "_q5_with_counts.rds"))
#         file_name_rds_raw <- file.path(output_dir, paste0("raw_results/act_ladel_fnorm_over_n", n, "_p", p, "_q5_raw_results.rds"))
#       } else if (data_gen_method == "q3noise") {
#         output_dir <- "~/multiDT/250409_simu/q3_noise"
#         dir.create(output_dir, showWarnings = FALSE)
#         dir.create(file.path(output_dir, "raw_results"), showWarnings = FALSE)
#         
#         file_name_excel <- file.path(output_dir, paste0("act_ladel_fnorm_noise_n", n, "_p", p, "_q3_with_counts.xlsx"))
#         file_name_rds_counts <- file.path(output_dir, paste0("act_ladel_fnorm_noise_n", n, "_p", p, "_q3_with_counts.rds"))
#         file_name_rds_raw <- file.path(output_dir, paste0("raw_results/act_ladel_fnorm_noise_n", n, "_p", p, "_q3_raw_results.rds"))
#       } else if (data_gen_method == "zip03") {
#         output_dir <- "~/multiDT/250409_simu/zip03"
#         dir.create(output_dir, showWarnings = FALSE)
#         dir.create(file.path(output_dir, "raw_results"), showWarnings = FALSE)
#         
#         file_name_excel <- file.path(output_dir, paste0("act_ladel_fnorm_zip03_n", n, "_p", p, "_q5_with_counts.xlsx"))
#         file_name_rds_counts <- file.path(output_dir, paste0("act_ladel_fnorm_zip03_n", n, "_p", p, "_q5_with_counts.rds"))
#         file_name_rds_raw <- file.path(output_dir, paste0("raw_results/act_ladel_fnorm_zip03_n", n, "_p", p, "_q5_raw_results.rds"))
#       } else if (data_gen_method == "unif") {
#         output_dir <- "~/multiDT/250409_simu/unif"
#         dir.create(output_dir, showWarnings = FALSE)
#         dir.create(file.path(output_dir, "raw_results"), showWarnings = FALSE)
#         
#         file_name_excel <- file.path(output_dir, paste0("act_ladel_funif_n", n, "_p", p, "_q5_with_counts.xlsx"))
#         file_name_rds_counts <- file.path(output_dir, paste0("act_ladel_funif_n", n, "_p", p, "_q5_with_counts.rds"))
#         file_name_rds_raw <- file.path(output_dir, paste0("raw_results/act_ladel_funif_n", n, "_p", p, "_q5_raw_results.rds"))
#       } else if (data_gen_method == "norm") {
#         output_dir <- "~/multiDT/250409_simu/norm"
#         dir.create(output_dir, showWarnings = FALSE)
#         dir.create(file.path(output_dir, "raw_results"), showWarnings = FALSE)
#         
#         file_name_excel <- file.path(output_dir, paste0("act_ladel_noiseN0-0.25_fnorm_n", n, "_p", p, "_q5_with_counts.xlsx"))
#         file_name_rds_counts <- file.path(output_dir, paste0("act_ladel_noiseN0-0.25_fnorm_n", n, "_p", p, "_q5_with_counts.rds"))
#         file_name_rds_raw <- file.path(output_dir, paste0("raw_results/act_ladel_noiseN0-0.25_fnorm_n", n, "_p", p, "_q5_raw_results.rds"))
#       } else if (data_gen_method == "q2noise") {
#         output_dir <- "~/multiDT/250409_simu/q2_noise"
#         dir.create(output_dir, showWarnings = FALSE)
#         dir.create(file.path(output_dir, "raw_results"), showWarnings = FALSE)
#         
#         file_name_excel <- file.path(output_dir, paste0("act_ladel_fnorm_noise_n", n, "_p", p, "_q2_with_counts.xlsx"))
#         file_name_rds_counts <- file.path(output_dir, paste0("act_ladel_fnorm_noise_n", n, "_p", p, "_q2_with_counts.rds"))
#         file_name_rds_raw <- file.path(output_dir, paste0("raw_results/act_ladel_fnorm_noise_n", n, "_p", p, "_q2_raw_results.rds"))
#       }
#       
#       # Create Excel workbook and save
#       wb <- createWorkbook()
#       addWorksheet(wb, "Results")
#       writeData(wb, "Results", final_data)
#       saveWorkbook(wb, file_name_excel, overwrite = TRUE)
#       
#       # Save count results as RDS
#       saveRDS(final_data, file_name_rds_counts)
#       
#       # Save raw results as RDS
#       saveRDS(all_results, file_name_rds_raw)
#       
#       cat("Files saved:\n")
#       cat("Excel:", file_name_excel, "\n")
#       cat("Counts RDS:", file_name_rds_counts, "\n")
#       cat("Raw Results RDS:", file_name_rds_raw, "\n")
#     }
#   }
#   
#   cat("\nSimulation for", data_gen_method, "completed.\n")
# }
# 
# # Example usage:
# # Run perbn simulation
# run_simulation(data_gen_method = "perbn", sample_sizes = c(100, 500), dimensions = c(100, 500), num_iterations = 10)
# 
# # Run other simulations (uncomment as needed)
#  run_simulation(data_gen_method = "perbp", sample_sizes = c(1000), dimensions = c(1000), num_iterations = 10)
# # run_simulation(data_gen_method = "over", sample_sizes = c(100, 500), dimensions = c(100, 500), num_iterations = 10)
# # run_simulation(data_gen_method = "q3noise", sample_sizes = c(100, 500), dimensions = c(100, 500), num_iterations = 10)
# # run_simulation(data_gen_method = "zip03", sample_sizes = c(100, 500), dimensions = c(100, 500), num_iterations = 10)
# # run_simulation(data_gen_method = "unif", sample_sizes = c(100, 500), dimensions = c(100, 500), num_iterations = 10)
# # run_simulation(data_gen_method = "norm", sample_sizes = c(100, 500), dimensions = c(100, 500), num_iterations = 10)
# # run_simulation(data_gen_method = "q2noise", sample_sizes = c(100, 500), dimensions = c(100, 500), num_iterations = 10)
# 
# # To run all simulations sequentially:
# # all_methods <- c("perbn", "perbp", "over", "q3noise", "zip03", "unif", "norm", "q2noise")
# # for (method in all_methods) {
# #   run_simulation(data_gen_method = method, sample_sizes = c(100, 500), dimensions = c(100, 500), num_iterations = 10)
# # }
#  
#  
#  
#  
#  
#  
# 
#  
#  S <- cov((log_centered_x)) 
#  eigvals_S <- eig(S)
#  R <- cor((log_centered_x), method = "pearson")
#  eigvals_R <- eig(R)
#  Rho <- cor((log_centered_x), method = "spearman")
#  eigvals_Rho <- eig(Rho)
#  MKen <- SK(t(log_centered_x))
#  eigvals_MKen <- eig(MKen)
#  ## Nadakuditi and Edelman (2008 tsp)
#  facnum_ne=NE(eigvals_S, p, n, maxfactor = 8)
#  ## Onatski (2010 restat)
#  facnum_ed=ED(eigvals_S, Kmax=8, 4)
#  ## Ahn and Horenstein (2013 eca)
#  facnum_gr=GR(eigvals_S, Kmax=8, p, n)
#  ## Lam and Yao (2012 aos)
#  facnum_er=ER(eigvals_MKen, Kmax=8)
#  
#  
 
 
 # Load required libraries
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
 library(ICtest)
 
 # Create cluster and register parallel backend if needed
 # cl <- makeCluster(20)
 # registerDoParallel(cl)
 # RhpcBLASctl::blas_set_num_threads(1)
 
 has_full_zero_row_or_col <- function(mat) {
   # Check if matrix has any full zero rows or columns
   any(rowSums(mat == 0) == ncol(mat)) || any(colSums(mat == 0) == nrow(mat))
 }
 
 # Function to run the simulation
 run_simulation <- function(data_gen_method = "perbn", sample_sizes = c(100, 500), dimensions = c(100, 500), num_iterations = 10) {
   # Create base output directory
   base_output_dir <- "~/multiDT/250409_simu"
   
   cat("Running simulation for", data_gen_method, "with s =", num_iterations, 
       ", n =", paste(sample_sizes, collapse=", "), 
       ", p =", paste(dimensions, collapse=", "), "\n")
   
   # Main Simulation Loop
   S <- num_iterations
   for (n in sample_sizes) {
     for (p in dimensions) {
       # Print current simulation parameters
       cat("\n=============================\n")
       cat("Scenario:", data_gen_method, "\n")
       cat("Sample size (n):", n, "\n")
       cat("Dimension (p):", p, "\n")
       cat("=============================\n")
       
       # Create storage for all results
       all_results <- list()
       # Initialize vectors with correct length
       facnum_ne_all <- numeric(S)
       facnum_ed_all <- numeric(S)
       facnum_gr_all <- numeric(S)
       facnum_er_all <- numeric(S)
       facnum_act_all <- numeric(S)
       facnum_ladel_all <- numeric(S)
       
       for (s in 1:S) {
         cat("Iteration", s, "of", S, "\n")
         
         # Source function file
         source("~/multiDT/250113function_tcv.R")
         
         # Set seed for reproducibility
         set.seed(s + 1000)
         
         # Factor matrix generation based on method
         if (data_gen_method == "q3noise" || data_gen_method == "q2noise") {
           # For q3noise and q2noise methods, set q accordingly
           q <- if(data_gen_method == "q3noise") 3 else 2
           FF <- matrix(rnorm(n * q), nrow = n, ncol = q)
         } else if (data_gen_method == "unif") {
           # For unif method, use uniform distribution for FF
           q <- 5
           FF <- matrix(runif(n * q, min = -2, max = 2), nrow = n, ncol = q)
         } else {
           # Default: use normal distribution with q=5
           q <- 5
           FF <- matrix(rnorm(n * q), nrow = n, ncol = q)
         }
         
         # Create a new list item for this iteration's results
         all_results[[s]] <- list()
         
         # Store q value in the current iteration's results
         all_results[[s]]$q <- q
         
         # Loading matrix generation (same for all methods)
         BB <- matrix(runif(p * q, min = -1, max = 1), nrow = p, ncol = q)
         
         # Additional parameter generation (same for all methods)
         a <- runif(p, min = 0, max = 1)
         
         # Identifiability adjustment
         result <- add_identifiability(FF, BB, a)
         FF0 <- result$H
         BB0 <- result$B
         alpha <- result$mu
         
         # Lambda parameter calculation based on method
         if (data_gen_method == "q3noise") {
           # Add noise to lambda for q3noise
           lambda <- exp(FF0 %*% t(BB0) + as.matrix(rep(1, n)) %*% t(as.matrix(alpha)) + 
                           matrix(rnorm(n * p, mean = 0, sd = 0.3), nrow = n, ncol = p))
         } else if (data_gen_method == "q2noise") {
           # Add noise to lambda for q2noise
           lambda <- exp(FF0 %*% t(BB0) + as.matrix(rep(1, n)) %*% t(as.matrix(alpha)) + 
                           matrix(rnorm(n * p, mean = 0, sd = 0.2), nrow = n, ncol = p))
         } else if (data_gen_method == "norm") {
           # Add noise to lambda for norm
           lambda <- exp(FF0 %*% t(BB0) + as.matrix(rep(1, n)) %*% t(as.matrix(alpha)) + 
                           matrix(rnorm(n * p, mean = 0, sd = 0.5), nrow = n, ncol = p))
         } else {
           # Default lambda calculation for other methods
           lambda <- exp(FF0 %*% t(BB0) + as.matrix(rep(1, n)) %*% t(as.matrix(alpha)))
         }
         
         # Store lambda in all_results
         all_results[[s]]$lambda <- lambda
         
         # Apply specific modifications based on method
         if (data_gen_method == "perbn") {
           # Perturb n (10% of the observations)
           lambda[((0.9*n+1):n),] <- 2*lambda[((0.9*n+1):n),]
         } else if (data_gen_method == "perbp") {
           # Perturb p (10% of the variables)
           lambda[,((0.9*p+1):p)] <- 2*lambda[,((0.9*p+1):p)]
         }
         
         # Data generation
         x <- matrix(nrow = n, ncol = p)
         
         if (data_gen_method == "over") {
           # Mixed distribution (overdispersion) data generation
           random_matrix <- runif(n*p)
           
           selection_threshold1 <- quantile(random_matrix, 0.1)
           selection_threshold2 <- quantile(random_matrix, 0.2)
           selection_threshold3 <- quantile(random_matrix, 0.3)
           
           for (i in 1:n) {
             for (j in 1:p) {
               if (random_matrix[(i-1)*p + j] <= selection_threshold1) {
                 # For random values <= threshold1, use theta = 0.1
                 theta <- 0.1
                 size <- theta
                 prob <- size / (size + lambda[i, j])
                 x[i, j] <- rnbinom(1, size, prob)
               } else if ((random_matrix[(i-1)*p + j] <= selection_threshold2) & 
                          (random_matrix[(i-1)*p + j] > selection_threshold1)) {
                 # For random values <= threshold2 but > threshold1, use theta = 1
                 theta <- 1
                 size <- theta
                 prob <- size / (size + lambda[i, j])
                 x[i, j] <- rnbinom(1, size, prob)
               } else if ((random_matrix[(i-1)*p + j] <= selection_threshold3) & 
                          (random_matrix[(i-1)*p + j] > selection_threshold2)) {
                 # For random values <= threshold3 but > threshold2, use theta = 10
                 theta <- 10
                 size <- theta
                 prob <- size / (size + lambda[i, j])
                 x[i, j] <- rnbinom(1, size, prob)
               } else {
                 # For random values > threshold3, generate Poisson distribution result
                 x[i, j] <- rpois(1, lambda[i, j])
               }
             }
           }
         } else {
           # Standard Poisson data generation for other methods
           for (i in 1:n) {
             for (j in 1:p) {
               x[i, j] <- rpois(1, lambda = lambda[i, j])
             }
           }
         }
         
         # Apply zero-inflation for zip03 method
         if (data_gen_method == "zip03") {
           repeat {
             # Create a copy of the original matrix to preserve non-zero structure
             x_temp <- x
             
             # Calculate number of elements to set to zero
             num_zeros <- ceiling(0.3 * length(x_temp))  # 30% of elements
             
             # Randomly select these element indices
             indices <- sample(length(x_temp), num_zeros, replace = FALSE)
             
             # Set selected elements to 0
             x_temp[indices] <- 0
             
             # Check if there are any full zero rows or columns
             if (!has_full_zero_row_or_col(x_temp)) {
               x <- x_temp  # Update x to new matrix
               break        # Exit loop
             }
           }
         }
         
         # Store generated data
         all_results[[s]]$FF <- FF
         all_results[[s]]$BB <- BB
         all_results[[s]]$a <- a
         all_results[[s]]$FF0 <- FF0
         all_results[[s]]$BB0 <- BB0
         all_results[[s]]$alpha <- alpha
         all_results[[s]]$x <- x
         
         # Preprocess data - centering log-transformed data
         log_centered_x <- sweep(log(x + 0.5), 1, rowMeans(log(x + 0.5)), "-")
         
         # Calculate matrices and eigenvalues for factor number estimation
         S_cov <- cov(log_centered_x)
         eigvals_S <- eig(S_cov)
         R <- cor(log_centered_x, method = "pearson")
         eigvals_R <- eig(R)
         Rho <- cor(log_centered_x, method = "spearman")
         eigvals_Rho <- eig(Rho)
         MKen <- SK(t(log_centered_x))
         eigvals_MKen <- eig(MKen)
         
         # 1. Nadakuditi and Edelman (2008 tsp)
         facnum_ne <- tryCatch({
           NE(eigvals_S, p, n, maxfactor = 8)
         }, error = function(e) {
           cat("Error occurred in NE (Nadakuditi and Edelman 2008):\n")
           print(e)
           cat("Current iteration:", s, "\n")
           NA
         })
         
         # 2. Onatski (2010 restat)
         facnum_ed <- tryCatch({
           ED(eigvals_S, Kmax = 8, 4)
         }, error = function(e) {
           cat("Error occurred in ED (Onatski 2010):\n")
           print(e)
           cat("Current iteration:", s, "\n")
           NA
         })
         
         # 3. Ahn and Horenstein (2013 eca)
         facnum_gr <- tryCatch({
           GR(eigvals_S, Kmax = 8, p, n)
         }, error = function(e) {
           cat("Error occurred in GR (Ahn and Horenstein 2013):\n")
           print(e)
           cat("Current iteration:", s, "\n")
           NA
         })
         
         # 4. Lam and Yao (2012 aos)
         facnum_er <- tryCatch({
           ER(eigvals_MKen, Kmax = 8,n,p)
         }, error = function(e) {
           cat("Error occurred in ER (Lam and Yao 2012):\n")
           print(e)
           cat("Current iteration:", s, "\n")
           NA
         })
         
         # 5. ACT method for factor number estimation
         facnum_act <- tryCatch({
           ACT(eigvals_R, (1 + sqrt(p / n)), Kmax = 8, n, p)
         }, error = function(e) {
           cat("Error occurred in ACT (Fan 2021 JASA):\n")
           print(e)
           cat("Current iteration:", s, "\n")
           NA
         })
         
         # 6. LADEL method for factor number estimation
         facnum_ladel <- tryCatch({
           #PCAladle(log_centered_x, ncomp = 8)$k
           1
         }, error = function(e) {
           cat("Error occurred in LADEL (Luo 2016):\n")
           print(e)
           cat("Current iteration:", s, "\n")
           NA
         })
         
         # Store factor number estimates
         facnum_ne_all[s] <- facnum_ne
         facnum_ed_all[s] <- facnum_ed
         facnum_gr_all[s] <- facnum_gr
         facnum_er_all[s] <- facnum_er
         facnum_act_all[s] <- facnum_act
         facnum_ladel_all[s] <- facnum_ladel
         
         all_results[[s]]$facnum_ne <- facnum_ne
         all_results[[s]]$facnum_ed <- facnum_ed
         all_results[[s]]$facnum_gr <- facnum_gr
         all_results[[s]]$facnum_er <- facnum_er
         all_results[[s]]$facnum_act <- facnum_act
         all_results[[s]]$facnum_ladel <- facnum_ladel
       }
       
       # Create results matrix with the new order
       facnum_mat <- cbind(
         facnum_ne = facnum_ne_all,
         facnum_ed = facnum_ed_all,
         facnum_gr = facnum_gr_all,
         facnum_er = facnum_er_all,
         facnum_act = facnum_act_all,
         facnum_ladel = facnum_ladel_all
       )
       
       # Convert to data frame
       res_df <- as.data.frame(facnum_mat)
       
       # Count occurrences of values 1 to 8 in each column
       count_data <- sapply(res_df, function(x) {
         tab <- table(factor(x, levels = 1:8))
         if (length(tab) < 8) {
           tab_full <- rep(0, 8)
           names(tab_full) <- as.character(1:8)
           tab_full[names(tab)] <- tab
           return(tab_full)
         } else {
           return(tab)
         }
       })
       
       # Combine the count data with the original results
       final_data <- rbind(count_data, res_df)
       
       # Prepare output directory and file names based on data generation method
       if (data_gen_method == "perbn") {
         output_dir <- "~/multiDT/250409_simu/perbn"
         dir.create(output_dir, showWarnings = FALSE)
         dir.create(file.path(output_dir, "raw_results"), showWarnings = FALSE)
         
         file_prefix <- paste0("add_methods_fnorm_perbn_n", n, "_p", p, "_q5")
         file_name_excel <- file.path(output_dir, paste0(file_prefix, "_with_counts.xlsx"))
         file_name_rds_counts <- file.path(output_dir, paste0(file_prefix, "_with_counts.rds"))
         file_name_rds_raw <- file.path(output_dir, "raw_results", paste0(file_prefix, "_raw_results.rds"))
       } else if (data_gen_method == "perbp") {
         output_dir <- "~/multiDT/250409_simu/perbp"
         dir.create(output_dir, showWarnings = FALSE)
         dir.create(file.path(output_dir, "raw_results"), showWarnings = FALSE)
         
         file_name_excel <- file.path(output_dir, paste0("add_methods_fnorm_perbp_n", n, "_p", p, "_q5_with_counts.xlsx"))
         file_name_rds_counts <- file.path(output_dir, paste0("add_methods_fnorm_perbp_n", n, "_p", p, "_q5_with_counts.rds"))
         file_name_rds_raw <- file.path(output_dir, paste0("raw_results/add_methods_fnorm_perbp_n", n, "_p", p, "_q5_raw_results.rds"))
       } else if (data_gen_method == "over") {
         output_dir <- "~/multiDT/250409_simu/over"
         dir.create(output_dir, showWarnings = FALSE)
         dir.create(file.path(output_dir, "raw_results"), showWarnings = FALSE)
         
         file_name_excel <- file.path(output_dir, paste0("add_methods_fnorm_over_n", n, "_p", p, "_q5_with_counts.xlsx"))
         file_name_rds_counts <- file.path(output_dir, paste0("add_methods_fnorm_over_n", n, "_p", p, "_q5_with_counts.rds"))
         file_name_rds_raw <- file.path(output_dir, paste0("raw_results/add_methods_fnorm_over_n", n, "_p", p, "_q5_raw_results.rds"))
       } else if (data_gen_method == "q3noise") {
         output_dir <- "~/multiDT/250409_simu/q3_noise"
         dir.create(output_dir, showWarnings = FALSE)
         dir.create(file.path(output_dir, "raw_results"), showWarnings = FALSE)
         
         file_name_excel <- file.path(output_dir, paste0("add_methods_fnorm_noise_n", n, "_p", p, "_q3_with_counts.xlsx"))
         file_name_rds_counts <- file.path(output_dir, paste0("add_methods_fnorm_noise_n", n, "_p", p, "_q3_with_counts.rds"))
         file_name_rds_raw <- file.path(output_dir, paste0("raw_results/add_methods_fnorm_noise_n", n, "_p", p, "_q3_raw_results.rds"))
       } else if (data_gen_method == "zip03") {
         output_dir <- "~/multiDT/250409_simu/zip03"
         dir.create(output_dir, showWarnings = FALSE)
         dir.create(file.path(output_dir, "raw_results"), showWarnings = FALSE)
         
         file_name_excel <- file.path(output_dir, paste0("add_methods_fnorm_zip03_n", n, "_p", p, "_q5_with_counts.xlsx"))
         file_name_rds_counts <- file.path(output_dir, paste0("add_methods_fnorm_zip03_n", n, "_p", p, "_q5_with_counts.rds"))
         file_name_rds_raw <- file.path(output_dir, paste0("raw_results/add_methods_fnorm_zip03_n", n, "_p", p, "_q5_raw_results.rds"))
       } else if (data_gen_method == "unif") {
         output_dir <- "~/multiDT/250409_simu/unif"
         dir.create(output_dir, showWarnings = FALSE)
         dir.create(file.path(output_dir, "raw_results"), showWarnings = FALSE)
         
         file_name_excel <- file.path(output_dir, paste0("add_methods_funif_n", n, "_p", p, "_q5_with_counts.xlsx"))
         file_name_rds_counts <- file.path(output_dir, paste0("add_methods_funif_n", n, "_p", p, "_q5_with_counts.rds"))
         file_name_rds_raw <- file.path(output_dir, paste0("raw_results/add_methods_funif_n", n, "_p", p, "_q5_raw_results.rds"))
       } else if (data_gen_method == "norm") {
         output_dir <- "~/multiDT/250409_simu/norm"
         dir.create(output_dir, showWarnings = FALSE)
         dir.create(file.path(output_dir, "raw_results"), showWarnings = FALSE)
         
         file_name_excel <- file.path(output_dir, paste0("add_methods_noiseN0-0.25_fnorm_n", n, "_p", p, "_q5_with_counts.xlsx"))
         file_name_rds_counts <- file.path(output_dir, paste0("add_methods_noiseN0-0.25_fnorm_n", n, "_p", p, "_q5_with_counts.rds"))
         file_name_rds_raw <- file.path(output_dir, paste0("raw_results/add_methods_noiseN0-0.25_fnorm_n", n, "_p", p, "_q5_raw_results.rds"))
       } else if (data_gen_method == "q2noise") {
         output_dir <- "~/multiDT/250409_simu/q2_noise"
         dir.create(output_dir, showWarnings = FALSE)
         dir.create(file.path(output_dir, "raw_results"), showWarnings = FALSE)
         
         file_name_excel <- file.path(output_dir, paste0("add_methods_fnorm_noise_n", n, "_p", p, "_q2_with_counts.xlsx"))
         file_name_rds_counts <- file.path(output_dir, paste0("add_methods_fnorm_noise_n", n, "_p", p, "_q2_with_counts.rds"))
         file_name_rds_raw <- file.path(output_dir, paste0("raw_results/add_methods_fnorm_noise_n", n, "_p", p, "_q2_raw_results.rds"))
       }
       
       # Create Excel workbook and save
       wb <- createWorkbook()
       addWorksheet(wb, "Results")
       writeData(wb, "Results", final_data)
       saveWorkbook(wb, file_name_excel, overwrite = TRUE)
       
       # Save count results as RDS
       saveRDS(final_data, file_name_rds_counts)
       print(final_data[1:8,])
       cat("\n")
       # Save raw results as RDS
       saveRDS(all_results, file_name_rds_raw)
       
       cat("Files saved:\n")
       cat("Excel:", file_name_excel, "\n")
       cat("Counts RDS:", file_name_rds_counts, "\n")
       cat("Raw Results RDS:", file_name_rds_raw, "\n")
     }
   }
   
   cat("\nSimulation for", data_gen_method, "completed.\n")
 }
 
cat("S=100, methods=all, n 100,200, 500, p 100,200, 500","\n")
 # To run all simulations sequentially:
 all_methods <- c("perbn", "norm", "over", "q3noise", "zip03", "unif", "perbp", "q2noise")
 for (method in all_methods) {
    run_simulation(data_gen_method = method, sample_sizes = c(100,200, 500), dimensions = c(100,200, 500), num_iterations = 100)
  }
 
 
 # Example usage:
 #run_simulation(data_gen_method = "perbn", sample_sizes = c(60,600), dimensions = c(60, 600), num_iterations = 10)
 # Run other simulations (uncomment as needed)
 #run_simulation(data_gen_method = "norm", sample_sizes = c(60,600), dimensions = c(60,600), num_iterations = 10)
 # run_simulation(data_gen_method = "over", sample_sizes = c(100, 500), dimensions = c(100, 500), num_iterations = 10)
 # run_simulation(data_gen_method = "q3noise", sample_sizes = c(100, 500), dimensions = c(100, 500), num_iterations = 10)
 #run_simulation(data_gen_method = "zip03", sample_sizes = c(60, 600), dimensions = c(60), num_iterations = 10)
 # run_simulation(data_gen_method = "unif", sample_sizes = c(100, 500), dimensions = c(100, 500), num_iterations = 10)
 # run_simulation(data_gen_method = "norm", sample_sizes = c(100, 500), dimensions = c(100, 500), num_iterations = 10)
 # run_simulation(data_gen_method = "q2noise", sample_sizes = c(100, 500), dimensions = c(100, 500), num_iterations = 10)
 
 
 