#### batch script for regressions on scf cluster ###

### SUBJECT IR16!!!! ###

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

## read in data ##
file_path_to_elecs_of_interest <- path(here(), "munge", "IR16_elecs_of_interest_insula.csv")
path_hg_clean <- path(here(), "munge", "IR16", "hg_behave_presentation_locked_extended_150.csv")
hg_behave <-  read.csv(path_hg_clean)
elecs_to_use <- read.csv(file_path_to_elecs_of_interest)

## regression parameters ##
# save info needed for regressions #
nBins <- colnames(hg_behave %>% select(starts_with("bin_")))
niter <- 10000

## electrodes ##
# prep #
all_elecs <- elecs_to_use %>% select(Electrode)
brain_behave_data <- hg_behave %>%
  filter(grepl(paste(all_elecs$Electrode, collapse = "|"), electrodes))
all_electrodes <- unique(brain_behave_data$electrodes)

## run regressions ##

# adv ineq #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_advent", nBins, region_name =  "Insula", niter, sub = "IR16", tag = "pres-locked")
# disadv ineq #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_disadvent", nBins, region_name =  "Insula", niter, sub = "IR16", tag = "pres-locked")
# self payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_payoff", nBins, region_name =  "Insula", niter, sub = "IR16", tag = "pres-locked")
# other payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_payoff", nBins, region_name =  "Insula", niter, sub = "IR16", tag = "pres-locked")
# self foregone #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_foregone", nBins, region_name =  "Insula", niter, sub = "IR16", tag = "pres-locked")
# other foregone #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_foregone", nBins, region_name =  "Insula", niter, sub = "IR16", tag = "pres-locked")
# self var paroff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_var_payoff", nBins, region_name =  "Insula", niter, sub = "IR16", tag = "pres-locked")
# other var payoff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_var_payoff", nBins, region_name =  "Insula", niter, sub = "IR16", tag = "pres-locked")
# self diff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_diff", nBins, region_name =  "Insula", niter, sub = "IR16", tag = "pres-locked")
# other diff #
run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_diff", nBins, region_name =  "Insula", niter, sub = "IR16", tag = "pres-locked")

