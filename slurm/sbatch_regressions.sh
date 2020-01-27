#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1

#SBATCH --mail-user=bstavel@berkeley.edu
#SBATCH --mail-type=ALL

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

R CMD BATCH --no-save ../analysis/batch_regression_script_par.R regression_ofc_1-17-19.out
