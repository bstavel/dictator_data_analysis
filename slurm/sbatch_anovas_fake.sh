#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=15
#SBATCH --nodes=1

#SBATCH --mail-user=bstavel@berkeley.edu
#SBATCH --mail-type=ALL

#export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

R CMD BATCH --no-save ../analysis/batch_stepped_anova_script_par.R stepped_anovas_fake_1-22-19.out
