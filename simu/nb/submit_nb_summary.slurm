#!/bin/bash -l
#SBATCH --job-name=TCV-NB-sum
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --output=logs/nb_summary_%j.out
#SBATCH --error=logs/nb_summary_%j.err

# Summary of all stages after the four jobs have finished (out_all_2510.R).
# 8 cores are requested only for memory: the script loads every checkpoint (about 6 GB for S = 100).
# Submit from this folder:  sbatch submit_nb_summary.slurm
IMG=/lustre/share/img/ood-app/rstudio_4.4.1.sif

cd "$SLURM_SUBMIT_DIR"
singularity exec $IMG Rscript --vanilla out_all_2510.R
