#!/bin/bash -l
#SBATCH --job-name=TCV-time
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --exclusive
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/time_%j.out
#SBATCH --error=logs/time_%j.err

# One job = one method script, on a node that no other job shares (--exclusive).
# --time: raise it if the partition allows; a job killed at the limit continues
# from its saved tasks when the same command is submitted again.
# Submit from the time/ directory after `mkdir -p logs`:
#   for m in TCV2 TCV5 TCV10 ICDT3 ICDT5 ICDT10; do sbatch submit_time.slurm $m; done
# TIME_WORKERS = number of concurrent single-threaded worker processes. The
# default 40 uses every core of the node (40 single-threaded workers); it must
# be the same for all six jobs. Override for all six with
#   TIME_WORKERS=20 sbatch submit_time.slurm TCV2   (and likewise for the others)
# If the partition mixes CPU models, add --constraint=<type> so that all six jobs
# run on the same node type; the hostname is recorded in every result row.
# After all six jobs have finished:  Rscript summarize_time.R
# Template adapted from submit_icdt.slurm (R 4.4.1 singularity image); verify the
# image path and the installed packages on the cluster before use.

set -e
METHOD=${1:?usage: sbatch submit_time.slurm TCV2|TCV5|TCV10|ICDT3|ICDT5|ICDT10}
R_IMAGE=/lustre/share/img/ood-app/rstudio_4.4.1.sif

# Single-threaded BLAS / OpenMP inside the container, set before R starts.
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1
export SINGULARITYENV_VECLIB_MAXIMUM_THREADS=1
export SINGULARITYENV_TIME_WORKERS=${TIME_WORKERS:-40}

cd "$SLURM_SUBMIT_DIR"
echo "$(date '+%F %T') START ${METHOD} on $(hostname) workers=${SINGULARITYENV_TIME_WORKERS} cpus=${SLURM_CPUS_PER_TASK}"
srun --ntasks=1 --cpus-per-task="${SLURM_CPUS_PER_TASK}" \
  singularity exec "$R_IMAGE" Rscript --vanilla "run_time_${METHOD}.R"
echo "$(date '+%F %T') FINISHED ${METHOD}"
