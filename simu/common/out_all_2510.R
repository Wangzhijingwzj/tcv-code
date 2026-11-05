
process_simulation_data <- function(folder_type, q_value = 5) {
  
  folder1 <- file.path("~/multiDT/250409_simu", folder_type)
  folder_save <- file.path("~/multiDT/250409_simu", folder_type, "output")
  
  
  dir.create(folder_save, recursive = TRUE, showWarnings = FALSE)
  
  
  n_values <- c(100, 200, 500, 1000)
  p_values <- c(100, 200, 500, 1000)
  
  
  orig_cols <- c("facnum_ratio", "facnum_bcv", "facnum_bcv2", "facnum_sr", 
                 "facnum_IC1", "facnum_IC2", "facnum_TCV2K", "facnum_TPCV2K", 
                 "facnum_TIPCV2K", "facnum_TCV5K", "facnum_TPCV5K", "facnum_TIPCV5K",
                 "facnum_TCV10K", "facnum_TPCV10K", "facnum_TIPCV10K")
  
  
  new_cols <- c("facnum_ne", "facnum_ed", "facnum_gr", "facnum_er", "facnum_ladel", "facnum_act")
  
  
  output_cols <- c("facnum_ne", "facnum_ed", "facnum_gr", "facnum_er", "facnum_ratio", 
                   "facnum_ladel", "facnum_act", "facnum_sr", "facnum_bcv", "facnum_bcv2",  
                   "facnum_IC1", "facnum_IC2", "facnum_TCV2K", "facnum_TPCV2K", 
                   "facnum_TIPCV2K", "facnum_TCV5K", "facnum_TPCV5K", "facnum_TIPCV5K",
                   "facnum_TCV10K", "facnum_TPCV10K", "facnum_TIPCV10K")
  
  
  file_patterns <- list(
    "norm" = list(
      orig_file = "noiseN0-0.25_fnorm_n%d_p%d_q%d_with_counts.xlsx",
      new_file = "add_methods_noiseN0-0.25_fnorm_n%d_p%d_q%d_with_counts.xlsx",
      output_file = "all_result2510_norm.xlsx",
      output_rds = "all_result2510_norm.rds"
    ),
    "over" = list(
      orig_file = "fnorm_over_n%d_p%d_q%d_with_counts.xlsx",
      new_file = "add_methods_fnorm_over_n%d_p%d_q%d_with_counts.xlsx",
      output_file = "all_result2510_over.xlsx",
      output_rds = "all_result2510_over.rds"
    ),
    "perbn" = list(
      orig_file = "fnorm_perbn_n%d_p%d_q%d_with_counts.xlsx",
      new_file = "add_methods_fnorm_perbn_n%d_p%d_q%d_with_counts.xlsx",
      output_file = "all_result2510_perbn.xlsx",
      output_rds = "all_result2510_perbn.rds"
    ),
    "perbp" = list(
      orig_file = "fnorm_perbp_n%d_p%d_q%d_with_counts.xlsx",
      new_file = "add_methods_fnorm_perbp_n%d_p%d_q%d_with_counts.xlsx",
      output_file = "all_result2510_perbp.xlsx",
      output_rds = "all_result2510_perbp.rds"
    ),
    "unif" = list(
      orig_file = "funif_n%d_p%d_q%d_with_counts.xlsx",
      new_file = "add_methods_funif_n%d_p%d_q%d_with_counts.xlsx",
      output_file = "all_result2510_unif.xlsx",
      output_rds = "all_result2510_unif.rds"
    ),
    "zip03" = list(
      orig_file = "fnorm_zip03_n%d_p%d_q%d_with_counts.xlsx",
      new_file = "add_methods_fnorm_zip03_n%d_p%d_q%d_with_counts.xlsx",
      output_file = "all_result2510_zip03.xlsx",
      output_rds = "all_result2510_zip03.rds"
    ),
    "q3_noise" = list(
      orig_file = "fnorm_noise_n%d_p%d_q%d_with_counts.xlsx",
      new_file = "add_methods_fnorm_noise_n%d_p%d_q%d_with_counts.xlsx",
      output_file = "all_result2510_q3_noise.xlsx",
      output_rds = "all_result2510_q3_noise.rds"
    ),
    "q2_noise" = list(
      orig_file = "fnorm_noise_n%d_p%d_q%d_with_counts.xlsx",
      new_file = "add_methods_fnorm_noise_n%d_p%d_q%d_with_counts.xlsx",
      output_file = "all_result2510_q2_noise.xlsx",
      output_rds = "all_result2510_q2_noise.rds"
    )
  )
  
  
  if (!folder_type %in% names(file_patterns)) {
    stop(sprintf("不支持的文件夹类型: %s", folder_type))
  }
  
  
  pattern <- file_patterns[[folder_type]]
  
  
  total_cols <- 2 + length(output_cols) + length(output_cols)
  result_matrix <- matrix(NA, nrow = length(n_values) * length(p_values), ncol = total_cols)
  
  
  colnames(result_matrix) <- c("n", "p", 
                               paste0(output_cols, "_correct(over|under)"),
                               paste0(output_cols, "_over_mean|under_mean")
  )
  
  row_index <- 1
  
  
  process_column <- function(dat, col, file_path, true_q) {
    
    S <- 100
    
    
    n_rows <- nrow(dat)
    
    
    if (n_rows >= (S+8)) {
      
      x <- dat[[col]][9:(S+8)]
    } else {
      
      x <- dat[[col]][9:n_rows]
      cat(sprintf("警告: %s 文件只有 %d 行，少于预期的 %d 行\n", 
                  basename(file_path), n_rows, (S+8)))
    }
    
    
    x <- x[!is.na(x)]
    
    x_nonzero <- x[x != 0]
    
    if (length(x) > length(x_nonzero)) {
      cat(sprintf("%s 文件列 %s 有 %d 个0元素，除去0元素一共有 %d 个元素\n", 
                  basename(file_path), col, length(x) - length(x_nonzero), length(x_nonzero)))
    }
    
    
    if (length(x_nonzero) > 0) {
      correct <- sum(floor(x_nonzero) == true_q, na.rm = TRUE)
      over <- sum(floor(x_nonzero) > true_q, na.rm = TRUE)
      under <- sum(floor(x_nonzero) < true_q, na.rm = TRUE)
      
      
      total_count <- correct + over + under
      if (total_count != length(x_nonzero)) {
        cat(sprintf("警告: %s 文件列 %s 的correct(%d)+over(%d)+under(%d)=%d，不等于非零元素个数%d\n", 
                    basename(file_path), col, correct, over, under, total_count, length(x_nonzero)))
      }
      
      
      over_mean <- 0
      if (over > 0) {
        over_values <- floor(x_nonzero)[floor(x_nonzero) > true_q]
        over_mean <- mean(over_values - true_q, na.rm = TRUE)
      }
      
      
      under_mean <- 0
      if (under > 0) {
        under_values <- floor(x_nonzero)[floor(x_nonzero) < true_q]
        under_mean <- mean(true_q - under_values, na.rm = TRUE)
      }
    } else {
      correct <- 0
      over <- 0
      under <- 0
      over_mean <- 0
      under_mean <- 0
    }
    
    correct_over_under <- sprintf("%d (%d|%d)", correct, over, under)
    over_under_mean <- sprintf("%.2f|%.2f", over_mean, under_mean)
    
    return(list(
      correct_over_under = correct_over_under,
      over_under_mean = over_under_mean
    ))
  }
  
  
  true_q <- q_value
  cat(sprintf("处理 %s 文件夹，真实因子个数为 %d\n", folder_type, true_q))
  
  
  for(n in n_values) {
    for(p in p_values) {
      
      file1 <- file.path(folder1, sprintf(pattern$orig_file, n, p, true_q))
      file2 <- file.path(folder1, sprintf(pattern$new_file, n, p, true_q))
      
      cat(sprintf("处理 n=%d, p=%d 的文件...\n", n, p))
      cat(sprintf("原始文件: %s\n", file1))
      cat(sprintf("新增文件: %s\n", file2))
      
      
      results <- rep(NA, total_cols)
      results[1:2] <- c(n, p)  
      
      
      file1_exists <- file.exists(file1)
      file2_exists <- file.exists(file2)
      
      if (!file1_exists) {
        cat(sprintf("警告: 原始文件不存在: %s\n", file1))
      }
      
      if (!file2_exists) {
        cat(sprintf("警告: 新增文件不存在: %s\n", file2))
      }
      
      
      if (!file1_exists && !file2_exists) {
        cat(sprintf("警告: n=%d, p=%d 的两个文件均不存在\n", n, p))
        slash_count <- total_cols - 2  
        result_matrix[row_index, ] <- c(n, p, rep("/", slash_count))
        row_index <- row_index + 1
        next
      }
      
      
      orig_results <- list()
      if (file1_exists) {
        tryCatch({
          cat(sprintf("正在读取原始文件: %s\n", basename(file1)))
          dat1 <- readxl::read_excel(file1)
          
          for(col in orig_cols) {
            if(col %in% colnames(dat1)) {
              cat(sprintf("处理列: %s\n", col))
              orig_results[[col]] <- process_column(dat1, col, file1, true_q)
            } else {
              cat(sprintf("警告: %s 文件中不存在列 %s\n", basename(file1), col))
              orig_results[[col]] <- list(
                correct_over_under = "NA (NA|NA)",
                over_under_mean = "NA|NA"
              )
            }
          }
        }, error = function(e) {
          cat(sprintf("错误: 处理 %s 文件时出错: %s\n", basename(file1), e$message))
          for(col in orig_cols) {
            orig_results[[col]] <- list(
              correct_over_under = "ERROR (ERROR|ERROR)",
              over_under_mean = "ERROR|ERROR"
            )
          }
        })
      } else {
        cat(sprintf("警告: n=%d, p=%d 的原始方法文件不存在\n", n, p))
        for(col in orig_cols) {
          orig_results[[col]] <- list(
            correct_over_under = "/ (/|/)",
            over_under_mean = "/|/"
          )
        }
      }
      
      
      new_results <- list()
      if (file2_exists) {
        tryCatch({
          cat(sprintf("正在读取新增文件: %s\n", basename(file2)))
          dat2 <- readxl::read_excel(file2)
          
          for(col in new_cols) {
            if(col %in% colnames(dat2)) {
              cat(sprintf("处理列: %s\n", col))
              new_results[[col]] <- process_column(dat2, col, file2, true_q)
            } else {
              cat(sprintf("警告: %s 文件中不存在列 %s\n", basename(file2), col))
              new_results[[col]] <- list(
                correct_over_under = "NA (NA|NA)",
                over_under_mean = "NA|NA"
              )
            }
          }
        }, error = function(e) {
          cat(sprintf("错误: 处理 %s 文件时出错: %s\n", basename(file2), e$message))
          for(col in new_cols) {
            new_results[[col]] <- list(
              correct_over_under = "ERROR (ERROR|ERROR)",
              over_under_mean = "ERROR|ERROR"
            )
          }
        })
      } else {
        cat(sprintf("警告: n=%d, p=%d 的新增方法文件不存在\n", n, p))
        for(col in new_cols) {
          new_results[[col]] <- list(
            correct_over_under = "/ (/|/)",
            over_under_mean = "/|/"
          )
        }
      }
      
      
      all_results <- c(orig_results, new_results)
      
      
      correct_over_under_vec <- c()
      over_under_mean_vec <- c()
      
      for(col in output_cols) {
        if(col %in% names(all_results)) {
          correct_over_under_vec <- c(correct_over_under_vec, all_results[[col]]$correct_over_under)
          over_under_mean_vec <- c(over_under_mean_vec, all_results[[col]]$over_under_mean)
        } else {
          cat(sprintf("警告: 结果中不存在列 %s\n", col))
          correct_over_under_vec <- c(correct_over_under_vec, "NA (NA|NA)")
          over_under_mean_vec <- c(over_under_mean_vec, "NA|NA")
        }
      }
      
      
      result_row <- c(n, p, correct_over_under_vec, over_under_mean_vec)
      result_matrix[row_index, ] <- result_row
      row_index <- row_index + 1
    }
  }
  
  
  result_df <- as.data.frame(result_matrix, stringsAsFactors = FALSE)
  
  
  colnames(result_df) <- c("n", "p", 
                           paste0(output_cols, "_correct(over|under)"),
                           paste0(output_cols, "_over_mean|under_mean"))
  
  
  xlsx_file <- file.path(folder_save, pattern$output_file)
  rds_file <- file.path(folder_save, pattern$output_rds)
  
  cat(sprintf("正在保存结果到: %s 和 %s\n", xlsx_file, rds_file))
  openxlsx::write.xlsx(result_df, file = xlsx_file)
  saveRDS(result_df, file = rds_file)
  
  cat(sprintf("%s 文件夹处理完成!\n\n", folder_type))
  
  return(list(
    xlsx_file = xlsx_file,
    rds_file = rds_file
  ))
}


required_packages <- c("readxl", "openxlsx")
for(pkg in required_packages) {
  if(!requireNamespace(pkg, quietly = TRUE)) {
    cat(sprintf("正在安装 %s 包...\n", pkg))
    install.packages(pkg)
  }
}


library(readxl)
library(openxlsx)


folder_types <- c(
  "norm", 
  "over", 
  "perbn", 
  "perbp", 
  "unif", 
  "zip03", 
  "q3_noise", 
  "q2_noise"
)


q_values <- c(
  "norm" = 5,
  "over" = 5,
  "perbn" = 5,
  "perbp" = 5,
  "unif" = 5,
  "zip03" = 5,
  "q3_noise" = 3,
  "q2_noise" = 2
)



results <- list()
for(folder in folder_types) {
  cat(sprintf("\n开始处理 %s 文件夹...\n", folder))
  results[[folder]] <- process_simulation_data(folder, q_values[folder])
}


cat("\n所有文件夹处理完成！结果如下：\n")
for(folder in folder_types) {
  cat(sprintf("%s: Excel文件位置 - %s\n", folder, results[[folder]]$xlsx_file))
  cat(sprintf("%s: RDS文件位置 - %s\n", folder, results[[folder]]$rds_file))
}
