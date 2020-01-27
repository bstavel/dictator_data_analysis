#### batch script for regressions on scf cluster ###

## libraries ##
library(tidyverse)
library(ggplot2)
library(magrittr)
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
source(path(here(), "R", "run_permuted_regressions_par.R"))
source(path(here(), "R", "compile_results.R"))
source(path(here(), "R", "step_wise_anova_method_par.R"))
source(path(here(), "R", "run_filtered_anova.R"))
source(path(here(), "R", 'mutate_cond.R'))

## paralellization ##
nCores <- 20
cl <- makeForkCluster(nCores)
registerDoParallel(cl)

## read in input data ##
path_hg_clean <- path(here(), "munge", "hg_behave.csv")
hg_behave <-  read.csv(path_hg_clean)

## read in regressions results data ##
regions_to_combine <- c("OFC", "Insula", "Cingulate", "STS")
regression_results <- compile_results(regions_to_combine)

## anovas ##
stepwise_anova_method_par(regression_results, hg_behave, niter = 1000)
