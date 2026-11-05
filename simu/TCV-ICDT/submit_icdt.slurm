#!/bin/bash -l
#SBATCH --job-name=TCV-ICDT
#SBATCH --partition=cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=40
#SBATCH --time=1-00:00:00
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=wangzhijing_wzj@163.com
#SBATCH --output=logs/icdt_%j.out
#SBATCH --error=logs/icdt_%j.err

set -e

# Submit from ~/multiDT/TCV-ICDT so that Slurm can open the logs/ files.
# Use the R 4.4.1 image from the user's existing server submission script.
ICDT_R_IMAGE=/lustre/share/img/ood-app/rstudio_4.4.1.sif

# Set these inside the container before R and its worker processes start.
export SINGULARITYENV_OMP_NUM_THREADS=1
export SINGULARITYENV_OPENBLAS_NUM_THREADS=1
export SINGULARITYENV_MKL_NUM_THREADS=1
export SINGULARITYENV_VECLIB_MAXIMUM_THREADS=1

cd "$HOME/multiDT/TCV-ICDT"

if [ "$#" -eq 0 ]; then
  set -- norm zip03 perbp over q3_noise unif perbn
fi

for scenario in "$@"; do
  echo "$(date '+%F %T') START ${scenario}"
  srun --ntasks=1 --cpus-per-task=40 \
    singularity exec "$ICDT_R_IMAGE" \
    Rscript --vanilla "${scenario}/run_icdt_${scenario}.R"
  echo "$(date '+%F %T') FINISHED ${scenario}"
done
