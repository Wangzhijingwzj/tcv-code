################################################################################
## 02_descriptive_h.R
## All numbers used in the reply to Reviewer #3, comment (h), on the real data.
## Input : data/peerj32_counts_44x130.rds   (made by 01_peerj32_matrix.R)
## Model : GFM::gfm(types = "poisson", q = 3, algorithm = "AM", maxIter = 30),
##         the three-factor fit of the paper
## Output: output/descriptive_h.txt, output/table_descriptive.tex,
##         output/fig_mean_variance.pdf, output/fig_mean_variance.png
## Sections
##   1 descriptive statistics: taxon means and variances, variance/mean, zeros, library sizes
##   2 where the zeros are (by taxon mean)
##   3 dispersion given the fitted intensities: taxon-wise negative binomial size,
##     profile likelihood with the intensities held fixed (descriptive only)
##   4 fit of the mean structure: proportion of Poisson deviance explained,
##     sum of squares, correlation of observed and fitted means
##   5 marginal variance implied by the fitted model (law of total variance):
##     Var(X_ij) = E(lambda_ij) + Var(lambda_ij), estimated by mean_i and var_i of the
##     fitted intensities, compared with the observed taxon variances
##   6 mean-variance figure
################################################################################
BASE <- "/Users/chengjunwang/Desktop/我的科研/TCV/csda/revision1/code/real-r1"; OUT <- file.path(BASE, "output")
dir.create(OUT, showWarnings = FALSE)
suppressPackageStartupMessages(library(GFM))
x <- readRDS(file.path(BASE, "data", "peerj32_counts_44x130.rds")); n <- nrow(x); p <- ncol(x)
mn <- colMeans(x); vr <- apply(x, 2, var); zr <- colMeans(x == 0); lib <- rowSums(x)
fit <- GFM::gfm(list(x), types = "poisson", q = 3, algorithm = "AM", verbose = FALSE)
lam <- exp(fit$hH %*% t(fit$hB) + rep(1, n) %*% t(fit$hmu))
q5 <- function(v) round(quantile(v, c(0, .25, .5, .75, 1)), 2)

sink(file.path(OUT, "descriptive_h.txt"), split = TRUE)
cat("Reviewer #3 (h): numbers used in the reply,", format(Sys.time()), "\n")
cat("data: peerj32, microbiome", as.character(packageVersion("microbiome")), "; n =", n, "samples, p =", p, "taxa; GFM", as.character(packageVersion("GFM")), "\n")

cat("\n== 1. descriptive statistics ==\n")
cat("library sizes:", range(lib), "\n")
cat("taxon mean      (min, Q1, median, Q3, max):", q5(mn), "\n")
cat("taxon variance  (min, Q1, median, Q3, max):", q5(vr), "\n")
cat("variance/mean   (min, Q1, median, Q3, max):", q5(vr / mn), "\n")
cat("proportion of zeros: overall", round(100 * mean(x == 0), 1), "% ; per taxon (min, Q1, median, Q3, max):", q5(zr), "\n")
sl <- summary(lm(log(vr) ~ log(mn)))$coefficients; cat(sprintf("slope of log variance on log mean: %.2f (se %.2f)\n", sl[2, 1], sl[2, 2]))

cat("\n== 2. where the zeros are ==\n")
z <- colSums(x == 0)
for (k in c(1, 10, 20)) cat(sprintf("taxa with mean < %-3d: %3d taxa, holding %.1f%% of all zeros\n", k, sum(mn < k), 100 * sum(z[mn < k]) / sum(z)))
cat(sprintf("taxa with mean >= 100: %d taxa, %d zeros among %d counts (%.2f%%)\n", sum(mn >= 100), sum(z[mn >= 100]), n * sum(mn >= 100), 100 * sum(z[mn >= 100]) / (n * sum(mn >= 100))))
cat(sprintf("for reference, a Poisson count with mean 0.5 / 1 / 5 is zero with probability %.1f%% / %.1f%% / %.1f%%\n", 100 * exp(-.5), 100 * exp(-1), 100 * exp(-5)))

cat("\n== 3. dispersion given the fitted intensities: taxon-wise NB size (profile likelihood, intensities fixed) ==\n")
profile_size <- function(y, mu) {
  z <- seq(-16, 20, length.out = 721); v <- vapply(z, function(s) sum(dnbinom(y, mu = mu, size = exp(s), log = TRUE)), numeric(1))
  i <- which.max(v); o <- optimize(function(s) -sum(dnbinom(y, mu = mu, size = exp(s), log = TRUE)), c(z[max(1, i - 1)], z[min(length(z), i + 1)]))
  if (sum(dpois(y, mu, log = TRUE)) >= -o$objective) Inf else exp(o$minimum)
}
sz <- sapply(1:p, function(j) profile_size(x[, j], lam[, j]))
cat("size (min, Q1, median, Q3, max):", q5(sz), "\n")
bins <- cut(sz, c(0, .1, 1, 10, Inf), right = FALSE, labels = c("< 0.1", "[0.1, 1)", "[1, 10)", ">= 10")); tb <- table(bins)
cat(paste(names(tb), sprintf("%d taxa (%.0f%%)", tb, 100 * tb / p), sep = ": ", collapse = " | "), "\n")
cat(sprintf("within [0.1, 10], the range of Example 3: %d taxa (%.0f%%)\n", sum(sz >= .1 & sz <= 10), 100 * mean(sz >= .1 & sz <= 10)))

cat("\n== 4. fit of the mean structure (in-sample) ==\n")
dev <- function(L) 2 * sum(ifelse(x > 0, x * log(x / L), 0) - (x - L)); m0 <- matrix(mn, n, p, byrow = TRUE)
cat(sprintf("Poisson deviance: taxon-means model %.0f, three-factor model %.0f, proportion explained %.1f%%\n", dev(m0), dev(lam), 100 * (1 - dev(lam) / dev(m0))))
cat(sprintf("sum of squares: %.3e -> %.3e, proportion explained %.1f%% | correlation of observed counts and fitted means %.3f\n", sum((x - m0)^2), sum((x - lam)^2), 100 * (1 - sum((x - lam)^2) / sum((x - m0)^2)), cor(as.vector(x), as.vector(lam))))

cat("\n== 5. marginal variance implied by the fitted model (law of total variance) ==\n")
ev <- colMeans(lam) + apply(lam, 2, var)
cat(sprintf("sum over taxa: observed %.3e, implied %.3e, ratio %.3f\n", sum(vr), sum(ev), sum(ev) / sum(vr)))
cat(sprintf("per-taxon ratio implied/observed: median %.2f; taxa with mean >= 1000: median %.2f\n", median(ev / vr), median((ev / vr)[mn >= 1000])))
cat(sprintf("implied variance/mean: median %.1f (observed %.1f)\n", median(ev / colMeans(lam)), median(vr / mn)))

cat("\n== 6. table and figure written ==\n")
fmt1 <- function(v) formatC(v, format = "f", digits = 1, big.mark = ",")
tex <- c("\\begin{table}[htbp]", "\\centering",
         sprintf("\\caption{Descriptive statistics of the gut microbiome data ($n=%d$ samples, $p=%d$ taxa). Quantiles are taken over the %d taxa. Library sizes range from %s to %s, and %.1f\\%% of all counts are zero. The last row is the taxon-wise negative binomial size estimated by profile likelihood with the intensities of the fitted three-factor model held fixed; it is a descriptive measure of extra-Poisson variation, not a refitted model.}",
                 n, p, p, format(min(lib), big.mark = ","), format(max(lib), big.mark = ","), 100 * mean(x == 0)),
         "\\label{tab-desc}", "\\begin{tabular}{lrrrrr}", "\\toprule", " & Min & Q1 & Median & Q3 & Max \\\\", "\\midrule",
         paste0("Taxon mean & ", paste(fmt1(quantile(mn, c(0, .25, .5, .75, 1))), collapse = " & "), " \\\\"),
         paste0("Taxon variance & ", paste(fmt1(quantile(vr, c(0, .25, .5, .75, 1))), collapse = " & "), " \\\\"),
         paste0("Variance/mean & ", paste(fmt1(quantile(vr / mn, c(0, .25, .5, .75, 1))), collapse = " & "), " \\\\"),
         paste0("Proportion of zeros & ", paste(formatC(quantile(zr, c(0, .25, .5, .75, 1)), format = "f", digits = 2), collapse = " & "), " \\\\"),
         paste0("NB size given the fitted means & ", paste(formatC(signif(quantile(sz, c(0, .25, .5, .75, 1)), 2), format = "g"), collapse = " & "), " \\\\"),
         "\\bottomrule", "\\end{tabular}", "\\end{table}")
writeLines(tex, file.path(OUT, "table_descriptive.tex"))
draw <- function() {
  blue <- "#2a78d6"; orange <- "#eb6834"
  par(mar = c(4.4, 5.4, 1, 1), las = 1, mgp = c(3.9, .7, 0))
  plot(mn, vr, log = "xy", type = "n", xaxt = "n", yaxt = "n", xlim = range(mn) * c(.7, 1.4), ylim = range(c(vr, ev, mn)) * c(.7, 1.4),
       xlab = "", ylab = "Taxon Variance")
  title(xlab = "Taxon Mean", line = 2.6)
  xt <- c(.1, 1, 10, 100, 1000); yt <- 10^(0:6)
  axis(1, at = xt, labels = c("0.1", "1", "10", "100", "1000"))
  axis(2, at = yt, labels = c("1", "10", "100", expression(10^3), expression(10^4), expression(10^5), expression(10^6)))
  abline(h = yt, v = xt, col = "#eeeeee"); abline(0, 1, lty = 2, col = "grey40")
  points(mn, vr, pch = 16, cex = .85, col = adjustcolor(blue, .85))
  points(mn, ev, pch = 17, cex = .85, col = adjustcolor(orange, .85))
  legend("topleft", bty = "n", pch = c(16, 17, NA), lty = c(NA, NA, 2), col = c(blue, orange, "grey40"),
         legend = c("Observed Value", "Fitted Value", "Variance = Mean"))
}
pdf(file.path(OUT, "fig_mean_variance.pdf"), width = 5.2, height = 4.6); draw(); dev.off()
png(file.path(OUT, "fig_mean_variance.png"), width = 5.2, height = 4.6, units = "in", res = 200, bg = "white"); draw(); dev.off()
cat("table_descriptive.tex, fig_mean_variance.pdf, fig_mean_variance.png\n")
sink()
