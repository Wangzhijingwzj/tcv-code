# README for tcv-code

# Determining the Number of Factors in Poisson Factor Models via Thinning Cross-Validation

This repository contains the scripts for the simulation studies, timing experiments, and real data analysis in the paper. The simulations and timing experiments were run on a Slurm cluster. For the real data analysis, we used `peerj32`, a dataset in the R package `microbiome` with 44 samples and 130 genus-level microbial taxa.

## Directory structure

**Shared functions**  
[main/](main/) contains the shared R functions for TCV and the C++ solver used in the real data analysis. [simu/common/](simu/common/) contains the C++ solver for the simulation studies, the script for the additional methods, and the results summary script.

**Simulation scenarios**  
The main scripts for Examples 1–4 are in [simu/norm/](simu/norm/), [simu/perbp/](simu/perbp/), [simu/over/](simu/over/), and [simu/zip03/](simu/zip03/), respectively. The main scripts for Supplementary Tables S2–S4 are in [simu/q3_noise/](simu/q3_noise/), [simu/unif/](simu/unif/), and [simu/perbn/](simu/perbn/), respectively.

**ICDT and negative binomial simulations**  
[simu/TCV-ICDT/](simu/TCV-ICDT/) contains the ICDT run scripts for the seven scenarios and the job submission script. [simu/nb/](simu/nb/) contains the negative binomial simulation for Table S5. It uses the wrapper and upstream implementation in [ICDT/](ICDT/) (`Wangzhijingwzj/ICDT`, commit `84aef20`).

**Timing and real data analysis**  
[time/](time/) contains the six timing scripts for Table 5 and the scripts that summarize the results and generate the table. [real/](real/) contains the scripts that prepare the data, select the number of factors in the Poisson model, and generate the figures.

**Paths**  
In the simulation scripts, `~/multiDT/250113function_tcv.R` refers to [main/250113function_tcv.R](main/250113function_tcv.R), and `~/glmcv/simu/confirm_jmle_omp_poisson_with_intercept_missing.cpp` refers to [the C++ file of the same name in simu/common/](simu/common/confirm_jmle_omp_poisson_with_intercept_missing.cpp). Create the output directory `~/multiDT/250409_simu/<scenario>/` before running each scenario. If you change the output root, use the same location in both the simulation and summary scripts.

**Slurm jobs**  
Before submitting ICDT, negative binomial, or timing jobs, run `mkdir -p logs` in the corresponding submission directory. Set the Singularity image path in each submission script to an R 4.4.1 image available on your cluster.

## Figures and tables in the main text

**Simulation grid**  
Each scenario uses all 16 combinations of `n,p ∈ {100,200,500,1000}`, with 100 repetitions per combination. As configured, `norm` runs `(1000,500)`, and `over` and `perbn` each run `(1000,1000)`. For `perbp`, add `(100,1000)`, `(200,1000)`, and `(500,1000)`. To cover the full grid, adjust the loops over n and p in these main scripts and extend `sample_sizes` and `dimensions` in the script for the additional methods to include all four values.

**Table 1 (Example 1)**  
Run [250410_TCV_norm.R](simu/norm/250410_TCV_norm.R). Then run the additional methods and summary scripts, in that order, to generate `norm/output/all_result2510_norm.xlsx`. With the saved simulation results in place, run [run_icdt_norm.R](simu/TCV-ICDT/norm/run_icdt_norm.R) to write the results for the ICDT column to `norm/results/summary_*/ICDT_summary.csv`.

**Table 2 (Example 2, feature contamination)**  
Run the batch scripts [250410_TCV_perbp.R](simu/perbp/250410_TCV_perbp.R), [250410_TCV_perbp1k125.R](simu/perbp/250410_TCV_perbp1k125.R), and [250410_TCV_perbp1k1k.R](simu/perbp/250410_TCV_perbp1k1k.R). Then run the additional methods and summary scripts, in that order, to generate `perbp/output/all_result2510_perbp.xlsx`. With the saved simulation results in place, run [run_icdt_perbp.R](simu/TCV-ICDT/perbp/run_icdt_perbp.R) to write the results for the ICDT column to `perbp/results/summary_*/ICDT_summary.csv`.

**Table 3 (Example 3, a mixture of Poisson and negative binomial counts)**  
Run [250410_TCV_over.R](simu/over/250410_TCV_over.R). Then run the additional methods and summary scripts, in that order, to generate `over/output/all_result2510_over.xlsx`. With the saved simulation results in place, run [run_icdt_over.R](simu/TCV-ICDT/over/run_icdt_over.R) to write the results for the ICDT column to `over/results/summary_*/ICDT_summary.csv`.

**Table 4 (Example 4, zero inflation)**  
Run the batch scripts [250410_TCV_zip03.R](simu/zip03/250410_TCV_zip03.R), [250410_TCV1251k_zip03.R](simu/zip03/250410_TCV1251k_zip03.R), [250410_TCV1k125_zip03.R](simu/zip03/250410_TCV1k125_zip03.R), and [250410_TCV1k1k_zip03.R](simu/zip03/250410_TCV1k1k_zip03.R). Then run the additional methods and summary scripts, in that order, to generate `zip03/output/all_result2510_zip03.xlsx`. With the saved simulation results in place, run [run_icdt_zip03.R](simu/TCV-ICDT/zip03/run_icdt_zip03.R) to write the results for the ICDT column to `zip03/results/summary_*/ICDT_summary.csv`.

**Table 5 (computing time)**  
From `time/`, use `sbatch submit_time.slurm <method>` to submit one job each for [TCV2](time/run_time_TCV2.R), [TCV5](time/run_time_TCV5.R), [TCV10](time/run_time_TCV10.R), [ICDT3](time/run_time_ICDT3.R), [ICDT5](time/run_time_ICDT5.R), and [ICDT10](time/run_time_ICDT10.R). When all six jobs have finished, run [summarize_time.R](time/summarize_time.R) and then [make_tables.R](time/make_tables.R) to generate `time/output/table_median_mad.tex` and `time/output/table_median_mad.csv`.

**Working directory for the real data analysis**  
Use `real/` as the working directory. In [V2real.R](real/V2real.R), point the absolute paths in `source()` and `sourceCpp()` to the R and C++ files of the same names in `main/`. The script reads `peerj32` directly and generates Figures 3 and 4, Table S6, and Figure S1 in a single run.

**Figure 2 (mean–variance plot) and descriptive statistics in Section 5.3**  
First run [01_peerj32_matrix.R](real/01_peerj32_matrix.R) to generate `real/data/peerj32_counts_44x130.rds`, then run [02_descriptive_h.R](real/02_descriptive_h.R). The outputs are `real/output/fig_mean_variance.pdf` and `real/output/descriptive_h.txt`.

**Figure 3 (scatterplots of factor loadings)**  
[V2real.R](real/V2real.R) writes `factor_loadings12_labeled_adjusted.pdf`, `factor_loadings13_labeled_adjusted.pdf`, and `factor_loadings23_labeled_adjusted.pdf` to the working directory.

**Figure 4 (biplots of sample scores and factor loadings)**  
[V2real.R](real/V2real.R) writes `biplot_F1_F2_keys_adjusted.pdf`, `biplot_F1_F3_keys_adjusted.pdf`, and `biplot_F2_F3_keys_adjusted.pdf` to the working directory.

## Supplementary figures and tables

**Table S1 (Example 1 with K = 2, 5, and 10)**  
Follow the sequence for Table 1. Run [250410_TCV_norm.R](simu/norm/250410_TCV_norm.R), then the script for the additional methods, and finally the summary script. Read the columns for the three values of K from `norm/output/all_result2510_norm.xlsx`.

**Table S2 (three factors in the data-generating model)**  
Run the batch scripts [250410_TCV_q3noise.R](simu/q3_noise/250410_TCV_q3noise.R), [250410_TCV1251k_q3noise.R](simu/q3_noise/250410_TCV1251k_q3noise.R), [250410_TCV1k125_q3noise.R](simu/q3_noise/250410_TCV1k125_q3noise.R), and [250410_TCV1k1k_q3noise.R](simu/q3_noise/250410_TCV1k1k_q3noise.R). Then run the additional methods and summary scripts, in that order, to generate `q3_noise/output/all_result2510_q3_noise.xlsx`. With the saved simulation results in place, run [run_icdt_q3_noise.R](simu/TCV-ICDT/q3_noise/run_icdt_q3_noise.R) to write the results for the ICDT column to `q3_noise/results/summary_*/ICDT_summary.csv`.

**Table S3 (factors drawn from the uniform distribution on (−2, 2))**  
Run the batch scripts [250410_TCV_unif.R](simu/unif/250410_TCV_unif.R), [250410_TCV1251k_unif.R](simu/unif/250410_TCV1251k_unif.R), [250410_TCV1k125_unif.R](simu/unif/250410_TCV1k125_unif.R), and [250410_TCV1k1k_unif.R](simu/unif/250410_TCV1k1k_unif.R). Then run the additional methods and summary scripts, in that order, to generate `unif/output/all_result2510_unif.xlsx`. With the saved simulation results in place, run [run_icdt_unif.R](simu/TCV-ICDT/unif/run_icdt_unif.R) to write the results for the ICDT column to `unif/results/summary_*/ICDT_summary.csv`.

**Table S4 (sample contamination)**  
Run [250410_TCV_perbn.R](simu/perbn/250410_TCV_perbn.R). Then run the additional methods and summary scripts, in that order, to generate `perbn/output/all_result2510_perbn.xlsx`. With the saved simulation results in place, run [run_icdt_perbn.R](simu/TCV-ICDT/perbn/run_icdt_perbn.R) to write the results for the ICDT column to `perbn/results/summary_*/ICDT_summary.csv`.

**Table S5 (negative binomial counts with size parameter 10)**  
In `simu/nb/`, run [250410_TCV_nb.R](simu/nb/250410_TCV_nb.R) and then [tcv_addmethods250506.R](simu/nb/tcv_addmethods250506.R). The ICDT script [run_icdt_nb.R](simu/nb/run_icdt_nb.R) can run in parallel with the main simulation. When all three scripts have finished, run [out_all_2510.R](simu/nb/out_all_2510.R) to generate `results_S100/output_eps1e20_*/comparison.csv` and `all_result2510_nb.xlsx` in the same output directory. In `comparison.csv`, the method names beginning with `facnum_NBest_` identify the three negative binomial criteria in Table S5.

**Submitting the jobs for Table S5**  
[nb_settings.R](simu/nb/nb_settings.R) splits the simulation grid into part A (p = 100, 200, 500) and part B (p = 1000). Submit [submit_nb_main.slurm](simu/nb/submit_nb_main.slurm) and [submit_nb_icdt.slurm](simu/nb/submit_nb_icdt.slurm) once with argument A and once with argument B. After all four jobs have finished, submit [submit_nb_summary.slurm](simu/nb/submit_nb_summary.slurm).

**Table S6 (criterion values for the real data)**  
[V2real.R](real/V2real.R) evaluates models with 1 to 10 factors for K = 2, 5, and 10, and writes `re_root_mat.xlsx` and `re_root_mat.csv`. The workbook groups the columns in the order K = 5, 2, 10. Use the column names to arrange the values for Table S6.

**Figure S1 (heatmap of the loading matrix)**  
[V2real.R](real/V2real.R) writes the heatmap of the loadings from the three-factor model to `factor_loadings_heatmap.pdf`.

## Dependencies

The original analyses used R 4.4.1 and GFM 1.2.1. Rcpp compiles the C++ source at run time. The solver in `main/` used for the real data analysis requires a compiler with OpenMP support.

**Model fitting, thinning, and comparison methods**  
`GFM`, `countsplit`, `bcv`, `ICtest`, `Matrix`, `irlba`, `pracma`, `pcaPP`, `Rfast`, `MASS`, `mvtnorm`, and `invgamma`.

**Parallel computing and the C++ interface**  
`parallel`, `foreach`, `doParallel`, `doSNOW`, `Rcpp`, and `RcppArmadillo`.

**Data processing and file I/O**  
`dplyr`, `tidyr`, `tibble`, `readr`, `readxl`, `openxlsx`, `stringr`, and `reshape2`. Using `getCurrentPath()` to locate a script in RStudio also requires `rstudioapi`.

**Real data analysis and plotting**  
`microbiome`, `phyloseq`, `curatedMetagenomicData`, `ggplot2`, `ggrepel`, `pheatmap`, `corrplot`, `RColorBrewer`, `grid`, `gridExtra`, `car`, `scatterplot3d`, and `ggpubr`.
