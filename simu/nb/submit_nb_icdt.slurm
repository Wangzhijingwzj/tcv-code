#!/bin/bash -l
#SBATCH --job-name=TCV-NB-icdt
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/nb_icdt_%j.out
#SBATCH --error=logs/nb_icdt_%j.err

# ICDT and VSDT on the same data (regenerated from the data seeds).
# Submit from this folder:  sbatch submit_nb_icdt.slurm A    (p = 100, 200, 500)
#                           sbatch submit_nb_icdt.slurm B    (p = 1000)
# n, p, S and the number of workers are set in nb_settings.R.
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1
IMG=/lustre/share/img/ood-app/rstudio_4.4.1.sif

cd "$SLURM_SUBMIT_DIR"
singularity exec $IMG Rscript --vanilla run_icdt_nb.R $1
