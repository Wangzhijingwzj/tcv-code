# ===============================================

# ===============================================
library(microbiome)
library(phyloseq)
library(dplyr)


data("peerj32")


cat("type:", class(peerj32), "\n")
print(peerj32)

data_peer=peerj32


sample_meta <- (data_peer$phyloseq@otu_table)
X.data=t(as.matrix(as.data.frame(sample_meta)))
dim(X.data)


Y=data_peer$meta$gender
length(Y)
#View(data_peer$meta)


library(curatedMetagenomicData)
library(dplyr)
library(ggplot2)


library(Rcpp)
library(RcppArmadillo)
library(parallel)
library(doParallel)
library(foreach)
library(dplyr)

Rcpp::sourceCpp("/Users/chengjunwang/Desktop/我的科研/TCV/confirm_jmle_omp_poisson_with_intercept_missing.cpp")

#Rcpp::sourceCpp("/Users/chengjunwang/Desktop/confirm_jmle_omp_poisson_with_intercept_missing.cpp")

#source("~/multiDT/250113function_tcv.R")
source("/Users/chengjunwang/Desktop/我的科研/TCV/250113function_tcv.R")
#X.data

library(foreach)
library(doParallel)
library(irlba)
library(doSNOW)
library(openxlsx)
library(countsplit)

library(readr)



n=dim(X.data)[1]
p=dim(X.data)[2]

re_TCV2=multiDT(X.data,K=2, rmax = 10)
log_TCV2K=re_TCV2$TCV
log_TCVPas2K=re_TCV2$TCV+seq(1:10)*((n+p)/(n*p))*log((n*p)/(n+p))
log_TICV2K=re_TCV2$TICV
log_TICVPas2K=re_TCV2$TICV+seq(1:10)*((n+p)/(n*p))*log((n*p)/(n+p))

## 3 factors
facnum_TCV2K=which.min(log_TCV2K)
cat("facnum_TCV2K = ",facnum_TCV2K,"\n")
cat("log_TCV2K = ",log_TCV2K,"\n")
facnum_TCVPas2K=which.min(log_TCVPas2K)
cat("facnum_TCVPas2K = ",facnum_TCVPas2K,"\n")
cat("log_TCVPas2K = ",log_TCVPas2K,"\n")
facnum_TICV2K=which.min(log_TICV2K)
cat("facnum_TICV2K = ",facnum_TICV2K,"\n")
cat("log_TICV2K = ",log_TICV2K,"\n")
facnum_TICVPas2K=which.min(log_TICVPas2K)
cat("facnum_TICVPas2K = ",facnum_TICVPas2K,"\n")
cat("log_TICVPas2K = ",log_TICVPas2K,"\n")


re_TCV5=multiDT(X.data,K=5, rmax = 10)
log_TCV5K=re_TCV5$TCV
log_TCVPas5K=re_TCV5$TCV+seq(1:10)*((n+p)/(n*p))*log((n*p)/(n+p))
log_TICV5K=re_TCV5$TICV
log_TICVPas5K=re_TCV5$TICV+seq(1:10)*((n+p)/(n*p))*log((n*p)/(n+p))

facnum_TCV5K=which.min(log_TCV5K)
facnum_TCVPas5K=which.min(log_TCVPas5K)
facnum_TICV5K=which.min(log_TICV5K)
facnum_TICVPas5K=which.min(log_TICVPas5K)


cat("facnum_TCV5K = ",facnum_TCV5K,"\n")
cat("log_TCV5K = ",log_TCV5K,"\n")
cat("facnum_TCVPas5K = ",facnum_TCVPas5K,"\n")
cat("log_TCVPas5K = ",log_TCVPas5K,"\n")
cat("facnum_TICV5K = ",facnum_TICV5K,"\n")
cat("log_TICV5K = ",log_TICV5K,"\n")
cat("facnum_TICVPas5K = ",facnum_TICVPas5K,"\n")
cat("log_TICVPas5K = ",log_TICVPas5K,"\n")





set.seed(1)
re_TCV10=multiDT(X.data,K=10, rmax = 10)
log_TCV10K=re_TCV10$TCV
log_TCVPas10K=re_TCV10$TCV+seq(1:10)*((n+p)/(n*p))*log((n*p)/(n+p))
log_TICV10K=re_TCV10$TICV
log_TICVPas10K=re_TCV10$TICV+seq(1:10)*((n+p)/(n*p))*log((n*p)/(n+p))

facnum_TCV10K=which.min(log_TCV10K)
facnum_TCVPas10K=which.min(log_TCVPas10K)
facnum_TICV10K=which.min(log_TICV10K)
facnum_TICVPas10K=which.min(log_TICVPas10K)

cat("facnum_TCV10K = ",facnum_TCV10K,"\n")
cat("log_TCV10K = ",log_TCV10K,"\n")
cat("facnum_TCVPas10K = ",facnum_TCVPas10K,"\n")
cat("log_TCVPas10K = ",log_TCVPas10K,"\n")
cat("facnum_TICV10K = ",facnum_TICV10K,"\n")
cat("log_TICV10K = ",log_TICV10K,"\n")
cat("facnum_TICVPas10K = ",facnum_TICVPas10K,"\n")
cat("log_TICVPas10K = ",log_TICVPas10K,"\n")



dim(X.data)
re_root_mat=as.matrix(cbind(log_TCV5K, log_TCVPas5K,log_TICVPas5K,log_TCV2K, log_TCVPas2K,log_TICVPas2K,log_TCV10K, log_TCVPas10K,log_TICVPas10K))
dim(re_root_mat)
print(re_root_mat)

write.csv(re_root_mat,file = "re_root_mat.csv")
library(openxlsx)
write.xlsx(re_root_mat,file = "re_root_mat.xlsx")



#############




re_3factors=(GFM::gfm(list(X.data), types = c("poisson"), maxIter = 30, q = 3, algorithm = "AM"))

temp_add=add_identifiability(re_3factors$hH,re_3factors$hB,re_3factors$hmu)
Fhat=temp_add$H
Lhat=temp_add$B
muhat=temp_add$mu


dim(Fhat)
dim(Lhat)
#View(Lhat)
rownames(Lhat)=colnames(X.data)
rownames(Fhat)=rownames(X.data)





library(ggplot2)
library(pheatmap)
library(corrplot)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(ggrepel)


print("Loading矩阵基本信息:")
print(paste("维度:", dim(Lhat)[1], "×", dim(Lhat)[2]))
print("前几个微生物的loading值:")
print(head(Lhat))

min(Lhat)
max(Lhat)

# > min(Lhat)
# [1] -2.875514
# > max(Lhat)
# [1] 4.223264


#pdf("factor_loadings_heatmap.pdf", width = 20, height = 10)
pheatmap(t(Lhat),  
         filename = "factor_loadings_heatmap.pdf", 
         cluster_rows = FALSE, 
         cluster_cols = TRUE,  
         color = colorRampPalette(c("blue", "white", "red"))(100),
         main = "Loadings Matrix Heatmap",
         fontsize_col = 9,  
         show_colnames = TRUE,  
         show_rownames = TRUE,  
         cellwidth = 8,   
         cellheight = 20,  
         legend_labels = NA)  
#dev.off()











# 
# library(ggplot2)
# library(reshape2)
# 






Lhat_long <- Lhat %>%
  as.data.frame() %>%
  mutate(Microbe = rownames(Lhat)) %>%
  pivot_longer(cols = -Microbe, names_to = "Factor", values_to = "Loading")





n_top <- 10  

for(i in 1:ncol(Lhat)) {
  factor_name <- paste("Factor", i)
  
  
  top_positive <- Lhat[order(Lhat[,i], decreasing = TRUE)[1:n_top], i, drop = FALSE]
  print(paste("\n", factor_name, "- 正向贡献最大的", n_top, "个微生物:"))
  print(top_positive)
  
  
  top_negative <- Lhat[order(Lhat[,i], decreasing = FALSE)[1:n_top], i, drop = FALSE]
  print(paste(factor_name, "- 负向贡献最大的", n_top, "个微生物:"))
  print(top_negative)
}







#total_loading <- sqrt(rowSums(Lhat^2))

total_loading=abs(Lhat[,1]) + abs(Lhat[,2])
extreme_microbes <- names(sort(total_loading, decreasing = TRUE)[1:10])

Lhat_df <- as.data.frame(Lhat)
Lhat_df$Microbe <- rownames(Lhat)
Lhat_df$is_extreme <- Lhat_df$Microbe %in% extreme_microbes





p6 <- ggplot(Lhat_df, aes(x = V1, y = V2)) +
  geom_point(aes(color = is_extreme, size = is_extreme), alpha = 0.6) +
  geom_text_repel(data = subset(Lhat_df, is_extreme), 
                  aes(label = Microbe), 
                  size = 6, 
                  max.overlaps = 15) +
  scale_color_manual(values = c("FALSE" = "lightgray", "TRUE" = "red")) +
  scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(x = "Factor 1 Loading", y = "Factor 2 Loading") +
  theme_minimal() +
  theme(
    legend.position = "none",
    
    axis.title = element_text(size = 18),
    
    axis.text = element_text(size = 16),
    
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )


ggsave("factor_loadings12_labeled_adjusted.pdf", p6, width = 10, height = 8)



###########################

total_loading=abs(Lhat[,1]) + abs(Lhat[,3])
extreme_microbes <- names(sort(total_loading, decreasing = TRUE)[1:10])

Lhat_df <- as.data.frame(Lhat)
Lhat_df$Microbe <- rownames(Lhat)
Lhat_df$is_extreme <- Lhat_df$Microbe %in% extreme_microbes

p7 <- ggplot(Lhat_df, aes(x = V1, y = V3)) +
  geom_point(aes(color = is_extreme, size = is_extreme), alpha = 0.6) +
  geom_text_repel(data = subset(Lhat_df, is_extreme), 
                  aes(label = Microbe),
                  size = 6, 
                  size = 3, max.overlaps = 15) +
  scale_color_manual(values = c("FALSE" = "lightgray", "TRUE" = "red")) +
  scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(x = "Factor 1 Loading", y = "Factor 3 Loading") +
  theme_minimal() +
  theme(
    legend.position = "none",
    
    axis.title = element_text(size = 18),
    
    axis.text = element_text(size = 16),
    
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )
ggsave("factor_loadings13_labeled_adjusted.pdf", p7, width = 10, height = 8)

######################


total_loading=abs(Lhat[,2]) + abs(Lhat[,3])
extreme_microbes <- names(sort(total_loading, decreasing = TRUE)[1:10])

Lhat_df <- as.data.frame(Lhat)
Lhat_df$Microbe <- rownames(Lhat)
Lhat_df$is_extreme <- Lhat_df$Microbe %in% extreme_microbes

p8 <- ggplot(Lhat_df, aes(x = V2, y = V3)) +
  geom_point(aes(color = is_extreme, size = is_extreme), alpha = 0.6) +
  geom_text_repel(data = subset(Lhat_df, is_extreme), 
                  aes(label = Microbe), 
                  size = 6, 
                  size = 3, max.overlaps = 15) +
  scale_color_manual(values = c("FALSE" = "lightgray", "TRUE" = "red")) +
  scale_size_manual(values = c("FALSE" = 1, "TRUE" = 2)) +
  geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.5) +
  geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.5) +
  labs(x = "Factor 2 Loading", y = "Factor 3 Loading") +
  theme_minimal() +
  theme(
    legend.position = "none",
    
    axis.title = element_text(size = 18),
    
    axis.text = element_text(size = 16),
    
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.5)
  )
ggsave("factor_loadings23_labeled_adjusted.pdf", p8, width = 10, height = 8)




##################
######


rownames(Fhat)=rownames(X.data)






library(ggplot2)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(gridExtra)
library(corrplot)
library(car)
#library(rgl)
library(scatterplot3d)
library(ggpubr)
library(reshape2)



factor_scores <- as.data.frame(Fhat)
colnames(factor_scores) <- c("Factor1", "Factor2", "Factor3")
factor_scores$Gender <- Y
factor_scores$Sample <- rownames(Fhat)






factor_long <- factor_scores %>%
  pivot_longer(cols = c(Factor1, Factor2, Factor3), 
               names_to = "Factor", values_to = "Score")



library(ggplot2)
library(ggrepel)
library(dplyr)
library(stringr)
library(grid)


factor_labels <- c(
  Factor1 = "Factor 1",
  Factor2 = "Factor 2",
  Factor3 = "Factor 3"
)
name=colnames(X.data)

loadings=as.data.frame(cbind(Lhat,name))
colnames(loadings)=c("Factor1","Factor2","Factor3","Microbe")
loadings <- loadings %>%
  mutate(across(c(Factor1, Factor2, Factor3), as.numeric))
gender_colors <- c("F"="#FF6B6B", "M"="#4ECDC4")


.select_keys_by_patterns <- function(loadings_df, patterns, ax_x, ax_y) {
  col_x <- paste0("Factor", ax_x)
  col_y <- paste0("Factor", ax_y)
  keep_list <- list()
  
  for (pat in patterns) {
    hit <- grepl(pat, loadings_df$Microbe, ignore.case = TRUE)
    if (!any(hit)) {
      warning(sprintf("Pattern '%s' not found in loadings$Microbe; skipped.", pat))
      next
    }
    sub <- loadings_df[hit, , drop = FALSE]
    mag <- sqrt(sub[[col_x]]^2 + sub[[col_y]]^2)
    pick <- sub[which.max(mag), , drop = FALSE]
    pick$Label <- pick$Microbe  
    keep_list[[length(keep_list) + 1]] <- pick
  }
  
  if (length(keep_list) == 0) stop("No key microbes matched any pattern.")
  bind_rows(keep_list) %>% distinct(Microbe, .keep_all = TRUE)
}


.compute_scale_factor <- function(load_sub, scores, ax_x, ax_y, pad = 0.85){
  rx <- max(abs(scores[[paste0("Factor", ax_x)]]), na.rm = TRUE)
  ry <- max(abs(scores[[paste0("Factor", ax_y)]]), na.rm = TRUE)
  lx <- max(abs(load_sub[[paste0("Factor", ax_x)]]), na.rm = TRUE)
  ly <- max(abs(load_sub[[paste0("Factor", ax_y)]]), na.rm = TRUE)
  min(pad * rx / ifelse(lx == 0, 1, lx),
      pad * ry / ifelse(ly == 0, 1, ly))
}



create_biplot_keys_patterns <- function(ax_x = 1, ax_y = 2,
                                        key_patterns,
                                        ellipse_level = 0.68,
                                        subtitle = NULL,
                                        gender_colors = c("F"="#FF6B6B", "M"="#4ECDC4")) {
  
  stopifnot(all(c("Microbe", paste0("Factor", ax_x), paste0("Factor", ax_y)) %in% names(loadings)))
  stopifnot(all(c(paste0("Factor", ax_x), paste0("Factor", ax_y), "Gender") %in% names(factor_scores)))
  
  
  load_sub <- .select_keys_by_patterns(loadings, key_patterns, ax_x, ax_y)
  
  
  sf <- .compute_scale_factor(load_sub, factor_scores, ax_x, ax_y, pad = 0.85)
  col_x <- paste0("Factor", ax_x); col_y <- paste0("Factor", ax_y)
  load_sub <- load_sub %>%
    mutate(xend = .data[[col_x]] * sf,
           yend = .data[[col_y]] * sf)
  
  
  x_lab <- factor_labels[[paste0("Factor", ax_x)]]
  if (is.null(x_lab)) x_lab <- paste("Factor", ax_x)
  y_lab <- factor_labels[[paste0("Factor", ax_y)]]
  if (is.null(y_lab)) y_lab <- paste("Factor", ax_y)
  
  
  p <- ggplot() +
    geom_point(data = factor_scores,
               aes(x = .data[[col_x]], y = .data[[col_y]], color = Gender),
               size = 3.2, alpha = 0.9) +
    stat_ellipse(data = factor_scores,
                 aes(x = .data[[col_x]], y = .data[[col_y]], color = Gender),
                 level = ellipse_level, linetype = 2, linewidth = 0.9, show.legend = FALSE) +
    geom_segment(data = load_sub,
                 aes(x = 0, y = 0, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.35, "cm"), type = "closed"),
                 linewidth = 1.15, color = "black", alpha = 0.9) +
    geom_text_repel(data = load_sub,
                    aes(x = xend, y = yend, label = Label),
                    size = 7, 
                    fontface = "bold", color = "black",
                    box.padding = 0.8, point.padding = 0.6, min.segment.length = 0) +
    geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.45) +
    geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.45) +
    scale_color_manual(values = gender_colors) +
    labs(subtitle = subtitle, x = x_lab, y = y_lab, color = "Gender") +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "top",
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey90"),
      plot.title = element_text(size = 16, face = "bold"),
      axis.title = element_text(size = 18), 
      axis.text = element_text(size = 16),  
      legend.title = element_text(size = 12),
      legend.text  = element_text(size = 11),
      
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.7)
    ) +
    coord_fixed()
  
  return(p)
}


keys_12 <- c("Mitsuokella", 
             "Prevotella melaninogenica", 
             "Eubacterium biforme",
             "^Campylobacter$", 
             "^Helicobacter$", 
             "Megamonas hypermegale")

p12 <- create_biplot_keys_patterns(ax_x = 1, ax_y = 2, key_patterns = keys_12)
ggsave("biplot_F1_F2_keys_adjusted.pdf", p12, width = 9, height = 8)






keys_13 <- c("Mitsuokella", 
             "Prevotella melaninogenica", 
             "Eubacterium biforme",
             "Escherichia coli", 
             "^Pseudomonas$", 
             "^Campylobacter$")

p13 <- create_biplot_keys_patterns(ax_x = 1, ax_y = 3, key_patterns = keys_13)
ggsave("biplot_F1_F3_keys_adjusted.pdf", p13, width = 9, height = 8)
#print(p13)




keys_23 <- c("Megamonas hypermegale", 
             "Lactobacillus catenaformis", 
             "Bacteroides vulgatus",
             "Prevotella melaninogenica", 
             "Escherichia coli", 
             "^Pseudomonas$")

p23 <- create_biplot_keys_patterns(ax_x = 2, ax_y = 3, key_patterns = keys_23)
ggsave("biplot_F2_F3_keys_adjusted.pdf", p23, width = 9, height = 8)
#print(p23)

