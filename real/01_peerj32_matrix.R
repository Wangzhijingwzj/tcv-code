################################################################################
## Run:  Rscript 01_peerj32_matrix.R
################################################################################
BASE <- "/Users/chengjunwang/Desktop/我的科研/TCV/csda/revision1/code/real-r1"
suppressPackageStartupMessages({ library(microbiome); library(phyloseq) })
data("peerj32")
data_peer <- peerj32
sample_meta <- (data_peer$phyloseq@otu_table)              # V2real.R line 18
X.data <- t(as.matrix(as.data.frame(sample_meta)))         # V2real.R line 19

dir.create(file.path(BASE, "data"), showWarnings = FALSE, recursive = TRUE)
write.csv(X.data, file.path(BASE, "data", "peerj32_counts_44x130.csv"))
saveRDS(X.data, file.path(BASE, "data", "peerj32_counts_44x130.rds"))
sig <- as.matrix(peerj32$microbes)[rownames(X.data), colnames(X.data)]
write.csv(sig, file.path(BASE, "data", "peerj32_hitchip_signal_44x130.csv"))
write.csv(microbiome::meta(peerj32$phyloseq), file.path(BASE, "data", "peerj32_sample_meta.csv"))

cat("matrix used in the paper: ", nrow(X.data), "samples x", ncol(X.data), "taxa\n")
cat("all integer:", all(X.data == round(X.data)), " min:", min(X.data), " max:", max(X.data),
    " zeros:", sum(X.data == 0), sprintf("(%.1f%%)", 100 * mean(X.data == 0)), "\n")
cat("packages: microbiome", as.character(packageVersion("microbiome")),
    ", phyloseq", as.character(packageVersion("phyloseq")), "\n")
cat("saved to", file.path(BASE, "data"), "\n")
cat("\nsamples (rows):", paste(rownames(X.data), collapse = " "), "\n")
cat("\ntaxa (columns):", paste(colnames(X.data), collapse = " | "), "\n")
cat("\npreview, first 6 samples x first 8 taxa:\n"); print(X.data[1:6, 1:8])
