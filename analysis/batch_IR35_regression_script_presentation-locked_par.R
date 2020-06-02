#### batch script for regressions on scf cluster ###

### SUBJECT IR35!!!! ###

## libraries ##
library(tidyverse)
library(ggplot2)
library(magrittr)
library(zoo)
library(lme4)
library(doParallel)
library(parallel)
library(foreach)
library(here)
library(fs)

## hand written functions ##
source(path(here(), "R", "load_behave_data.R"))
source(path(here(), "R", "prep_behave_data.R"))
source(path(here(), "R", "stretch_start_end.R"))
source(path(here(), "R", "load_high_gamma_data.R"))
source(path(here(), "R", "rolling_window_and_baseline.R"))
source(path(here(), "R", "run_permuted_regressions.R"))
source(path(here(), "R", "run_permuted_regressions_par.R"))
source(path(here(), "R", "run_filtered_anova.R"))
source(path(here(), "R", 'mutate_cond.R'))

## paralellization ##
nCores <- 10
registerDoParallel(nCores)

## regression parameters ##
# save info needed for regressions #
niter <- 10000

## electrodes ##
# load #
file_path_to_elecs_of_interest <- path(here(), "munge", "IR35_elecs_of_interest_insula.csv")
elecs_to_use <- read.csv(file_path_to_elecs_of_interest)
# prep #
all_elecs <- elecs_to_use %>% select(Electrode)

### HFA ###

## read in data ##
path_hg_clean <- path(here(), "munge", "IR35", "hfa_behave_presentation_locked_wavelets.csv")
hfa_behave <-  read.csv(path_hg_clean)
# merge with elecs #
brain_behave_data <- hfa_behave %>%
  filter(grepl(paste(all_elecs$Electrode, collapse = "|"), electrodes))
all_electrodes <- unique(brain_behave_data$electrodes)
# bin names #
nBins <- colnames(hfa_behave %>% select(starts_with("bin_")))

## run regressions ##

# adv ineq #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_advent", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# disadv ineq #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_disadvent", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# self payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_payoff", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# other payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_payoff", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# self foregone #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_foregone", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# other foregone #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_foregone", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# self var paroff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_var_payoff", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# other var payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_var_payoff", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# # self diff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_diff", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# # other diff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_diff", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")

### theta ###

## read in data ##
path_hg_clean <- path(here(), "munge", "IR35", "theta_behave_presentation_locked_wavelets.csv")
theta_behave <-  read.csv(path_hg_clean)
# merge with elecs #
brain_behave_data <- theta_behave %>%
  filter(grepl(paste(all_elecs$Electrode, collapse = "|"), electrodes))
all_electrodes <- unique(brain_behave_data$electrodes)
# bin names #
nBins <- colnames(theta_behave %>% select(starts_with("bin_")))

## run regressions ##

# adv ineq #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_advent", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# disadv ineq #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_disadvent", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# self payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_payoff", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# other payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_payoff", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# self foregone #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_foregone", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# other foregone #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_foregone", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# self var paroff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_var_payoff", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# other var payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_var_payoff", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# # self diff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_diff", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")
# # other diff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_diff", nBins, region_name =  "Insula", niter, sub = "IR35", tag = "pres-locked-wavelets")

