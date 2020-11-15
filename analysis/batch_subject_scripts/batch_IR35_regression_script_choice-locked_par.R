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
nCores <- 32
registerDoParallel(nCores)

## regression parameters ##
# save info needed for regressions #
niter <- 1000

## electrodes ##
# load #
file_path_to_elecs_of_interest <- path(here(), "munge", "IR35_elecs_of_interest.csv")
elecs_to_use <- read.csv(file_path_to_elecs_of_interest)
# prep #
all_elecs <- elecs_to_use %>% select(Electrode)

### theta ###

# ## read in data ##
# path_hp_clean <- path(here(), "munge", "IR35", "theta_behave_choice_rscaler_2575_200.csv")
# power_behave <-  read.csv(path_hp_clean)
# # merge with elecs #
# brain_behave_data <- power_behave %>%
#   filter(grepl(paste(all_elecs$Electrode, collapse = "|"), electrodes))
# all_electrodes <- unique(brain_behave_data$electrodes)
# # bin names #
# nBins <- colnames(brain_behave_data %>% select(starts_with("pre_"), starts_with("post_")))

## run regressions ##

# # adv ineq #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_advent", nBins, region_name =  "All", niter, sub = "IR35", tag = "theta-choice-locked-hilbertRS")
# # disadv ineq #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_disadvent", nBins, region_name =  "All", niter, sub = "IR35", tag = "theta-choice-locked-hilbertRS")
# # self payoff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "theta-choice-locked-hilbertRS")
# # other payoff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "theta-choice-locked-hilbertRS")
# # self foregone #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_foregone", nBins, region_name =  "All", niter, sub = "IR35", tag = "theta-choice-locked-hilbertRS")
# # other foregone #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_foregone", nBins, region_name =  "All", niter, sub = "IR35", tag = "theta-choice-locked-hilbertRS")
# # self var paroff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_var_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "theta-choice-locked-hilbertRS")
# # other var payoff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_var_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "theta-choice-locked-hilbertRS")
# # self diff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_diff", nBins, region_name =  "All", niter, sub = "IR35", tag = "theta-choice-locked-hilbertRS")
# # other diff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_diff", nBins, region_name =  "All", niter, sub = "IR35", tag = "theta-choice-locked-hilbertRS")


### HFA ###

# ## read in data ##
# path_hp_clean <- path(here(), "munge", "IR35", "hfa_behave_choice_rscaler_2575_200.csv")
# power_behave <-  read.csv(path_hp_clean)
# # merge with elecs #
# brain_behave_data <- power_behave %>%
#   filter(grepl(paste(all_elecs$Electrode, collapse = "|"), electrodes))
# all_electrodes <- unique(brain_behave_data$electrodes)
# # bin names #
# nBins <- colnames(brain_behave_data %>% select(starts_with("pre_"), starts_with("post_")))

## run regressions ##

# # adv ineq #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_advent", nBins, region_name =  "All", niter, sub = "IR35", tag = "hfa-choice-locked-hilbertRS")
# # disadv ineq #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_disadvent", nBins, region_name =  "All", niter, sub = "IR35", tag = "hfa-choice-locked-hilbertRS")
# # self payoff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "hfa-choice-locked-hilbertRS")
# # other payoff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "hfa-choice-locked-hilbertRS")
# # self foregone #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_foregone", nBins, region_name =  "All", niter, sub = "IR35", tag = "hfa-choice-locked-hilbertRS")
# # other foregone #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_foregone", nBins, region_name =  "All", niter, sub = "IR35", tag = "hfa-choice-locked-hilbertRS")
# # self var paroff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_var_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "hfa-choice-locked-hilbertRS")
# # other var payoff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_var_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "hfa-choice-locked-hilbertRS")
# # self diff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_diff", nBins, region_name =  "All", niter, sub = "IR35", tag = "hfa-choice-locked-hilbertRS")
# # other diff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_diff", nBins, region_name =  "All", niter, sub = "IR35", tag = "hfa-choice-locked-hilbertRS")


### beta ###

## read in data ##
path_hp_clean <- path(here(), "munge", "IR35", "beta_behave_choice_rscaler_2575_200.csv")
power_behave <-  read.csv(path_hp_clean)
# merge with elecs #
brain_behave_data <- power_behave %>%
  filter(grepl(paste(all_elecs$Electrode, collapse = "|"), electrodes))
all_electrodes <- unique(brain_behave_data$electrodes)
# bin names #
nBins <- colnames(brain_behave_data %>% select(starts_with("pre_"), starts_with("post_")))

## run regressions ##

# # adv ineq #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_advent", nBins, region_name =  "All", niter, sub = "IR35", tag = "beta-choice-locked-hilbertRS")
# disadv ineq #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_disadvent", nBins, region_name =  "All", niter, sub = "IR35", tag = "beta-choice-locked-hilbertRS")
# self payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "beta-choice-locked-hilbertRS")
# other payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "beta-choice-locked-hilbertRS")
# # self foregone #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_foregone", nBins, region_name =  "All", niter, sub = "IR35", tag = "beta-choice-locked-hilbertRS")
# # other foregone #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_foregone", nBins, region_name =  "All", niter, sub = "IR35", tag = "beta-choice-locked-hilbertRS")
# # self var paroff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_var_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "beta-choice-locked-hilbertRS")
# # other var payoff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_var_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "beta-choice-locked-hilbertRS")
# # self diff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_diff", nBins, region_name =  "All", niter, sub = "IR35", tag = "beta-choice-locked-hilbertRS")
# # other diff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_diff", nBins, region_name =  "All", niter, sub = "IR35", tag = "beta-choice-locked-hilbertRS")


### gamma ###

## read in data ##
path_hp_clean <- path(here(), "munge", "IR35", "gamma_behave_choice_rscaler_2575_200.csv")
power_behave <-  read.csv(path_hp_clean)
# merge with elecs #
brain_behave_data <- power_behave %>%
  filter(grepl(paste(all_elecs$Electrode, collapse = "|"), electrodes))
all_electrodes <- unique(brain_behave_data$electrodes)
# bin names #
nBins <- colnames(brain_behave_data %>% select(starts_with("pre_"), starts_with("post_")))

## run regressions ##

# adv ineq #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_advent", nBins, region_name =  "All", niter, sub = "IR35", tag = "gamma-choice-locked-hilbertRS")
# disadv ineq #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_disadvent", nBins, region_name =  "All", niter, sub = "IR35", tag = "gamma-choice-locked-hilbertRS")
# self payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "gamma-choice-locked-hilbertRS")
# other payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "gamma-choice-locked-hilbertRS")
# # self foregone #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_foregone", nBins, region_name =  "All", niter, sub = "IR35", tag = "gamma-choice-locked-hilbertRS")
# # other foregone #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_foregone", nBins, region_name =  "All", niter, sub = "IR35", tag = "gamma-choice-locked-hilbertRS")
# # self var paroff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_var_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "gamma-choice-locked-hilbertRS")
# # other var payoff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_var_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "gamma-choice-locked-hilbertRS")
# # self diff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_diff", nBins, region_name =  "All", niter, sub = "IR35", tag = "gamma-choice-locked-hilbertRS")
# # other diff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_diff", nBins, region_name =  "All", niter, sub = "IR35", tag = "gamma-choice-locked-hilbertRS")

### alpha ###

## read in data ##
path_hp_clean <- path(here(), "munge", "IR35", "alpha_behave_choice_rscaler_2575_200.csv")
power_behave <-  read.csv(path_hp_clean)
# merge with elecs #
brain_behave_data <- power_behave %>%
  filter(grepl(paste(all_elecs$Electrode, collapse = "|"), electrodes))
all_electrodes <- unique(brain_behave_data$electrodes)
# bin names #
nBins <- colnames(brain_behave_data %>% select(starts_with("pre_"), starts_with("post_")))

## run regressions ##

# adv ineq #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_advent", nBins, region_name =  "All", niter, sub = "IR35", tag = "alpha-choice-locked-hilbertRS")
# disadv ineq #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_disadvent", nBins, region_name =  "All", niter, sub = "IR35", tag = "alpha-choice-locked-hilbertRS")
# self payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "alpha-choice-locked-hilbertRS")
# other payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "alpha-choice-locked-hilbertRS")
# # self foregone #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_foregone", nBins, region_name =  "All", niter, sub = "IR35", tag = "alpha-choice-locked-hilbertRS")
# # other foregone #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_foregone", nBins, region_name =  "All", niter, sub = "IR35", tag = "alpha-choice-locked-hilbertRS")
# # self var paroff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_var_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "alpha-choice-locked-hilbertRS")
# # other var payoff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_var_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "alpha-choice-locked-hilbertRS")
# # self diff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_diff", nBins, region_name =  "All", niter, sub = "IR35", tag = "alpha-choice-locked-hilbertRS")
# # other diff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_diff", nBins, region_name =  "All", niter, sub = "IR35", tag = "alpha-choice-locked-hilbertRS")

### delta ###

## read in data ##
path_hp_clean <- path(here(), "munge", "IR35", "delta_behave_choice_rscaler_2575_200.csv")
power_behave <-  read.csv(path_hp_clean)
# merge with elecs #
brain_behave_data <- power_behave %>%
  filter(grepl(paste(all_elecs$Electrode, collapse = "|"), electrodes))
all_electrodes <- unique(brain_behave_data$electrodes)
# bin names #
nBins <- colnames(brain_behave_data %>% select(starts_with("pre_"), starts_with("post_")))

## run regressions ##

# adv ineq #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_advent", nBins, region_name =  "All", niter, sub = "IR35", tag = "delta-choice-locked-hilbertRS")
# disadv ineq #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_disadvent", nBins, region_name =  "All", niter, sub = "IR35", tag = "delta-choice-locked-hilbertRS")
# self payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "delta-choice-locked-hilbertRS")
# other payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "delta-choice-locked-hilbertRS")
# # self foregone #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_foregone", nBins, region_name =  "All", niter, sub = "IR35", tag = "delta-choice-locked-hilbertRS")
# # other foregone #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_foregone", nBins, region_name =  "All", niter, sub = "IR35", tag = "delta-choice-locked-hilbertRS")
# # self var paroff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_var_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "delta-choice-locked-hilbertRS")
# # other var payoff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_var_payoff", nBins, region_name =  "All", niter, sub = "IR35", tag = "delta-choice-locked-hilbertRS")
# # self diff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_diff", nBins, region_name =  "All", niter, sub = "IR35", tag = "delta-choice-locked-hilbertRS")
# # other diff #
# run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_diff", nBins, region_name =  "All", niter, sub = "IR35", tag = "delta-choice-locked-hilbertRS")
