#### batch script for regressions on scf cluster ###

## libraries ##
library(tidyverse)
library(ggplot2)
library(magrittr)
# library(ggthemr)
library(grid)
library(gtable)
library(gridExtra)
library(wesanderson)
library(ggsci)
library(zoo)
library(kableExtra)
library(lme4)
library(RColorBrewer)
library(doParallel)
library(parallel)
library(foreach)
library(here)
library(fs)
library(ggcorrplot)

## hand written functions ##
source(path(here(), "R", "load_behave_data.R"))
source(path(here(), "R", "prep_behave_data.R"))
source(path(here(), "R", "stretch_start_end.R"))
source(path(here(), "R", "load_high_gamma_data.R"))
source(path(here(), "R", "rolling_window_and_baseline.R"))
source(path(here(), "R", "run_permuted_regressions.R"))
source(path(here(), "R", "run_permuted_regressions_par.R"))
source(path(here(), "R", "run_permuted_regressions.R"))
source(path(here(), "R", "run_filtered_anova.R"))
# source('~/Projects/nice_r_functions/ggpaired_pretty.R')
source(path(here(), "R", 'mutate_cond.R'))

## plotting helpers ##
# ggthemr("solarized")

## paralellization ##
nCores <- 10
registerDoParallel(nCores)


## read in data ##
file_path_to_elecs_of_interest <- path(here(), "munge", "IR35_elecs_of_interest.csv")
path_hg_clean <- path(here(), "munge", "hg_behave.csv")
hg_behave <-  read.csv(path_hg_clean)
elecs_to_use <- read.csv(file_path_to_elecs_of_interest)

## regressions ##
# save info needed for regressions #
nBins <- colnames(hg_behave %>% select(starts_with("bin_")))
niter <- 10000 

## ofc ##
# prep #
ofc_elecs <- elecs_to_use %>% filter(ofc_loc_meeting == 1) %>% select(Electrode)
brain_behave_data <- hg_behave %>% 
  filter(grepl(paste(ofc_elecs$Electrode, collapse = "|"), electrodes))
ofc_electrodes <- unique(brain_behave_data$electrodes)

## run regressions ##
# adv ineq #
run_permuted_regressions_par(brain_behave_data, electrodes = ofc_electrodes, regressor = "ineq_advent", nBins, region_name = "OFC", niter)
# disadv ineq #
# run_permuted_regressions(brain_behave_data, electrodes = ofc_electrodes, nBins, region_name = "OFC", niter)

# run anovas #
# run_filtered_anova(brain_behave_data, region_name = "OFC")

# ## insula ##
# # prep #
# insula_elecs <- elecs_to_use %>% filter(insula_loc_meeting == 1) %>% select(Electrode)
# brain_behave_data <- hg_behave %>% 
#   filter(grepl(paste(insula_elecs$Electrode, collapse = "|"), electrodes))
# insula_electrodes <- unique(brain_behave_data$electrodes)
# # regresions #
# run_permuted_regressions(brain_behave_data, electrodes = insula_electrodes, nBins, region_name = "Insula", niter)
# # anova #
# run_filtered_anova(brain_behave_data, region_name = "Insula")
# 
# 
# # # mfg #
# # mfg_elecs <- elecs_to_use %>% filter(mfg_loc_meeting == 1) %>% select(Electrode)
# # brain_behave_data <- hg_behave %>% 
# #   filter(grepl(paste(mfg_elecs$Electrode, collapse = "|"), electrodes))
# # mfg_electrodes <- unique(brain_behave_data$electrodes)
# # 
# # run_permuted_regressions(brain_behave_data, electrodes = mfg_electrodes, nBins, region_name = "MFG", niter)
# 
# ## sts ##
# # prep #
# sts_elecs <- elecs_to_use %>% filter(sts_loc_meeting == 1) %>% select(Electrode)
# brain_behave_data <- hg_behave %>% 
#   filter(grepl(paste(sts_elecs$Electrode, collapse = "|"), electrodes))
# sts_electrodes <- unique(brain_behave_data$electrodes)
# # regressions #
# run_permuted_regressions(brain_behave_data, electrodes = sts_electrodes, nBins, region_name = "STS", niter)
# # anova #
# run_filtered_anova(brain_behave_data, region_name = "STS")
# 
