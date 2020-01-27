#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=20
#SBATCH --nodes=1

#SBATCH --mail-user=bstavel@berkeley.edu
#SBATCH --mail-type=ALL


R CMD BATCH --no-save ../analysis/batch_stepped_anova_script_par.R stepped_anovas_1-24-19.out
