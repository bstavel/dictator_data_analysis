#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1

#SBATCH --mail-user=bstavel@berkeley.edu
#SBATCH --mail-type=ALL

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_all_subs_multiple_regression_script_pres-locked_par.R  multiple_regressions_presentation_all_12-30-2020.out

R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_all_subs_multiple_regression_script_choice-locked_par.R  multiple_regressions_choice_all_12-30-2020.out

#R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_IR28_regression_script_presentation-locked_par.R  regressions_presentation_IR28_10-19-2020.out

#R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_IR28_regression_script_choice-locked_par.R  regressions_choice_IR28_10-19-2020.out

#R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_IR35_regression_script_presentation-locked_par.R  regressions_presentation_IR35_12-21-2020.out

#R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_IR35_regression_script_choice-locked_par.R  regressions_choice_IR35_12-21-2020.out

#R CMD BATCH --no-save ../analysis/batch_IR19_regression_script_presentation-locked_par.R  regressions_pres_IR19_7-22-2020.out

#R CMD BATCH --no-save ../analysis/batch_IR19_regression_script_choice-locked_par.R  regressions_choice_IR19_7-28-2020.out

#R CMD BATCH --no-save ../analysis/batch_IR39_regression_script_presentation-locked_par.R  regressions_pres_IR39_4-26-19.out

#R CMD BATCH --no-save ../analysis/batch_IR39_regression_script_choice-locked_par.R  regressions_choice_IR39_6-3-2020.out

#R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_IR57_regression_script_presentation-locked_par.R  regressions_pres_IR57_10-19-2020.out

#R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_IR57_regression_script_choice-locked_par.R  regressions_choice_IR57_10-19-2020.out

#R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_IR9_regression_script_presentation-locked_par.R  regressions_pres_IR9_11-11-2020.out

#R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_IR9_regression_script_choice-locked_par.R  regressions_choice_IR9_11-1-2020.out

#R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_IR10_regression_script_presentation-locked_par.R  regressions_pres_IR10_11-1-2020.out

#R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_IR10_regression_script_choice-locked_par.R  regressions_choice_IR10_11-1-2020.out

#R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_IR26_regression_script_presentation-locked_par.R  regressions_pres_IR26_11-1-2020.out

#R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_IR26_regression_script_choice-locked_par.R  regressions_choice_IR26_11-1-2020.out

#R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_IR16_regression_script_presentation-locked_par.R  regressions_pres_IR16_11-1-2020.out

#R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_IR16_regression_script_choice-locked_par.R  regressions_choice_IR16_11-1-2020.out

#R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_CP34_regression_script_presentation-locked_par.R  regressions_pres_CP34_11-23-2020.out

#R CMD BATCH --no-save ../analysis/batch_subject_scripts/batch_CP34_regression_script_choice-locked_par.R  regressions_choice_CP34_11-23-2020.out
