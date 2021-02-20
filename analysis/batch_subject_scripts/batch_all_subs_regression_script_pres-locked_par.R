#### batch script for regressions on scf cluster ###

### ALL SUBJECTS!!!! ###

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
source(path(here(), "R", "run_permuted_regressions_par.R"))
source(path(here(), "R", "run_filtered_anova.R"))
source(path(here(), "R", 'mutate_cond.R'))

## paralellization ##
nCores <- 20
registerDoParallel(nCores)

## regression parameters ##
# save info needed for regressions #
niter <- 1000

## subs to run ##
subs <- c("IR9", "IR10", "IR16", "IR26", "IR28", "IR35", "IR57", "CP34") 

for(sub in subs){
  # Null main file out each time #
  power_behave <- NULL

    ## electrodes ##
    # load #
    file_path_to_elecs_of_interest <- path(here(), "munge", paste0(sub, "_elecs_of_interest.csv"))
    elecs_to_use <- read.csv(file_path_to_elecs_of_interest)
    # prep #
    all_elecs <- elecs_to_use %>% select(Electrode)
    
    ### theta ###
    
    ## read in data ##
    path_hp_clean <- path(here(), "munge", sub, "theta_behave_presentation_rscaler_2575_200.csv")
    power_behave <-  read.csv(path_hp_clean)
    # merge with elecs #
    brain_behave_data <- power_behave %>%
      filter(grepl(paste(all_elecs$Electrode, collapse = "|"), electrodes))
    all_electrodes <- unique(brain_behave_data$electrodes)
    # bin names #
    nBins <- colnames(power_behave %>% select(starts_with("bin_")))
    
    ## run regressions ##
    # adv ineq #
    run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_var", nBins, region_name =  "All", niter, sub = sub, tag = "theta-pres-locked-hilbertRS")
    # # disadv ineq #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_disadvent", nBins, region_name =  "All", niter, sub = sub, tag = "theta-pres-locked-hilbertRS")
    # # self payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "theta-pres-locked-hilbertRS")
    # # other payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "theta-pres-locked-hilbertRS")
    # # self foregone #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_foregone", nBins, region_name =  "All", niter, sub = sub, tag = "theta-pres-locked-hilbertRS")
    # # other foregone #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_foregone", nBins, region_name =  "All", niter, sub = sub, tag = "theta-pres-locked-hilbertRS")
    # # self var paroff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_var_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "theta-pres-locked-hilbertRS")
    # # other var payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_var_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "theta-pres-locked-hilbertRS")
    # # self diff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_diff", nBins, region_name =  "All", niter, sub = sub, tag = "theta-pres-locked-hilbertRS")
    # # other diff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_diff", nBins, region_name =  "All", niter, sub = sub, tag = "theta-pres-locked-hilbertRS")
    
    
    ### HFA ###

    ## read in data ##
    path_hp_clean <- path(here(), "munge", sub, "hfa_behave_presentation_rscaler_2575_200.csv")
    power_behave <-  read.csv(path_hp_clean)
    # merge with elecs #
    brain_behave_data <- power_behave %>%
      filter(grepl(paste(all_elecs$Electrode, collapse = "|"), electrodes))
    all_electrodes <- unique(brain_behave_data$electrodes)
    # bin names #
    nBins <- colnames(power_behave %>% select(starts_with("bin_")))

    ## run regressions ##

    # adv ineq #
    run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_var", nBins, region_name =  "All", niter, sub = sub, tag = "hfa-pres-locked-hilbertRS")
    # # disadv ineq #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_disadvent", nBins, region_name =  "All", niter, sub = sub, tag = "hfa-pres-locked-hilbertRS")
    # # self payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "hfa-pres-locked-hilbertRS")
    # # other payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "hfa-pres-locked-hilbertRS")
    # # self foregone #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_foregone", nBins, region_name =  "All", niter, sub = sub, tag = "hfa-pres-locked-hilbertRS")
    # # other foregone #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_foregone", nBins, region_name =  "All", niter, sub = sub, tag = "hfa-pres-locked-hilbertRS")
    # # self var paroff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_var_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "hfa-pres-locked-hilbertRS")
    # # other var payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_var_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "hfa-pres-locked-hilbertRS")
    # # # self diff #
    # # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_diff", nBins, region_name =  "All", niter, sub = sub, tag = "hfa-pres-locked-hilbertRS")
    # # # other diff #
    # # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_diff", nBins, region_name =  "All", niter, sub = sub, tag = "hfa-pres-locked-hilbertRS")


    ### beta ###

    ## read in data ##
    path_hp_clean <- path(here(), "munge", sub, "beta_behave_presentation_rscaler_2575_200.csv")
    power_behave <-  read.csv(path_hp_clean)
    # merge with elecs #
    brain_behave_data <- power_behave %>%
      filter(grepl(paste(all_elecs$Electrode, collapse = "|"), electrodes))
    all_electrodes <- unique(brain_behave_data$electrodes)
    # bin names #
    nBins <- colnames(power_behave %>% select(starts_with("bin_")))

    ## run regressions ##

    # # adv ineq #
    run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_var", nBins, region_name =  "All", niter, sub = sub, tag = "beta-pres-locked-hilbertRS")
    # # disadv ineq #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_disadvent", nBins, region_name =  "All", niter, sub = sub, tag = "beta-pres-locked-hilbertRS")
    # # self payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "beta-pres-locked-hilbertRS")
    # # other payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "beta-pres-locked-hilbertRS")
    # # self foregone #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_foregone", nBins, region_name =  "All", niter, sub = sub, tag = "beta-pres-locked-hilbertRS")
    # # other foregone #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_foregone", nBins, region_name =  "All", niter, sub = sub, tag = "beta-pres-locked-hilbertRS")
    # # self var paroff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_var_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "beta-pres-locked-hilbertRS")
    # # other var payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_var_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "beta-pres-locked-hilbertRS")
    # # self diff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_diff", nBins, region_name =  "All", niter, sub = sub, tag = "beta-pres-locked-hilbertRS")
    # # other diff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_diff", nBins, region_name =  "All", niter, sub = sub, tag = "beta-pres-locked-hilbertRS")


    ### gamma ###

    ## read in data ##
    path_hp_clean <- path(here(), "munge", sub, "gamma_behave_presentation_rscaler_2575_200.csv")
    power_behave <-  read.csv(path_hp_clean)
    # merge with elecs #
    brain_behave_data <- power_behave %>%
      filter(grepl(paste(all_elecs$Electrode, collapse = "|"), electrodes))
    all_electrodes <- unique(brain_behave_data$electrodes)
    # bin names #
    nBins <- colnames(power_behave %>% select(starts_with("bin_")))

    ## run regressions ##

    # adv ineq #
    run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_var", nBins, region_name =  "All", niter, sub = sub, tag = "gamma-pres-locked-hilbertRS")
    # # disadv ineq #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_disadvent", nBins, region_name =  "All", niter, sub = sub, tag = "gamma-pres-locked-hilbertRS")
    # # self payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "gamma-pres-locked-hilbertRS")
    # # other payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "gamma-pres-locked-hilbertRS")
    # # self foregone #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_foregone", nBins, region_name =  "All", niter, sub = sub, tag = "gamma-pres-locked-hilbertRS")
    # # other foregone #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_foregone", nBins, region_name =  "All", niter, sub = sub, tag = "gamma-pres-locked-hilbertRS")
    # # self var paroff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_var_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "gamma-pres-locked-hilbertRS")
    # # other var payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_var_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "gamma-pres-locked-hilbertRS")
    # # self diff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_diff", nBins, region_name =  "All", niter, sub = sub, tag = "gamma-pres-locked-hilbertRS")
    # # other diff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_diff", nBins, region_name =  "All", niter, sub = sub, tag = "gamma-pres-locked-hilbertRS")

    ### alpha ###

    ## read in data ##
    path_hp_clean <- path(here(), "munge", sub, "alpha_behave_presentation_rscaler_2575_200.csv")
    power_behave <-  read.csv(path_hp_clean)
    # merge with elecs #
    brain_behave_data <- power_behave %>%
      filter(grepl(paste(all_elecs$Electrode, collapse = "|"), electrodes))
    all_electrodes <- unique(brain_behave_data$electrodes)
    # bin names #
    nBins <- colnames(power_behave %>% select(starts_with("bin_")))

    ## run regressions ##

    # adv ineq #
    run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_var", nBins, region_name =  "All", niter, sub = sub, tag = "alpha-pres-locked-hilbertRS")
    # # disadv ineq #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_disadvent", nBins, region_name =  "All", niter, sub = sub, tag = "alpha-pres-locked-hilbertRS")
    # # self payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "alpha-pres-locked-hilbertRS")
    # # other payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "alpha-pres-locked-hilbertRS")
    # # self foregone #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_foregone", nBins, region_name =  "All", niter, sub = sub, tag = "alpha-pres-locked-hilbertRS")
    # # other foregone #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_foregone", nBins, region_name =  "All", niter, sub = sub, tag = "alpha-pres-locked-hilbertRS")
    # # self var paroff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_var_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "alpha-pres-locked-hilbertRS")
    # # other var payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_var_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "alpha-pres-locked-hilbertRS")
    # # self diff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_diff", nBins, region_name =  "All", niter, sub = sub, tag = "alpha-pres-locked-hilbertRS")
    # # other diff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_diff", nBins, region_name =  "All", niter, sub = sub, tag = "alpha-pres-locked-hilbertRS")

    ### delta ###

    ## read in data ##
    path_hp_clean <- path(here(), "munge", sub, "delta_behave_presentation_rscaler_2575_200.csv")
    power_behave <-  read.csv(path_hp_clean)
    # merge with elecs #
    brain_behave_data <- power_behave %>%
      filter(grepl(paste(all_elecs$Electrode, collapse = "|"), electrodes))
    all_electrodes <- unique(brain_behave_data$electrodes)
    # bin names #
    nBins <- colnames(power_behave %>% select(starts_with("bin_")))

    ## run regressions ##

    # adv ineq #
    run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_var", nBins, region_name =  "All", niter, sub = sub, tag = "delta-pres-locked-hilbertRS")
    # # disadv ineq #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "ineq_disadvent", nBins, region_name =  "All", niter, sub = sub, tag = "delta-pres-locked-hilbertRS")
    # # self payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "delta-pres-locked-hilbertRS")
    # # other payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "delta-pres-locked-hilbertRS")
    # # self foregone #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_foregone", nBins, region_name =  "All", niter, sub = sub, tag = "delta-pres-locked-hilbertRS")
    # # other foregone #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_foregone", nBins, region_name =  "All", niter, sub = sub, tag = "delta-pres-locked-hilbertRS")
    # # self var paroff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_var_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "delta-pres-locked-hilbertRS")
    # # other var payoff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_var_payoff", nBins, region_name =  "All", niter, sub = sub, tag = "delta-pres-locked-hilbertRS")
    # # self diff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "self_diff", nBins, region_name =  "All", niter, sub = sub, tag = "delta-pres-locked-hilbertRS")
    # # other diff #
    # run_permuted_regressions_par(brain_behave_data, electrodes = all_electrodes, regressor = "other_diff", nBins, region_name =  "All", niter, sub = sub, tag = "delta-pres-locked-hilbertRS")

}