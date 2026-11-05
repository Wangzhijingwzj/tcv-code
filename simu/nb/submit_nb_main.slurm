#!/bin/bash -l
#SBATCH --job-name=TCV-NB-main
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --time=2-00:00:00
#SBATCH --output=logs/nb_main_%j.out
#SBATCH --error=logs/nb_main_%j.err

# Main simulation, then the additional methods (SE, ED, GR, ER, ACT).
# Submit from this folder:  sbatch submit_nb_main.slurm A    (p = 100, 200, 500)
#                           sbatch submit_nb_main.slurm B    (p = 1000)
# n, p, S and the number of workers are set in nb_settings.R.
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1
IMG=/lustre/share/img/ood-app/rstudio_4.4.1.sif

cd "$SLURM_SUBMIT_DIR"
singularity exec $IMG Rscript --vanilla 250410_TCV_nb.R $1 &&
singularity exec $IMG Rscript --vanilla tcv_addmethods250506.R $1
