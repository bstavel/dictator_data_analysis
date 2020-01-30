#### batch script for regressions on scf cluster ###

## libraries ##
if(!require('tidyverse')){
  install.packages("tidyverse", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}
if(!require('ggplot2')){
  install.packages("ggplot2", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}
if(!require('magrittr')){
  install.packages("magrittr", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}
if(!require('grid')){
  install.packages("grid", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}
if(!require('gtable')){
  install.packages("gtable", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}
if(!require('gridExtra')){
  install.packages("gridExtra", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}
if(!require('wesanderson')){
  install.packages("wesanderson", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}
if(!require('zoo')){
  install.packages("zoo", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}
if(!require('kableExtra')){
  install.packages("kableExtra", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}
if(!require('lme4')){
  install.packages("lme4", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}
if(!require('RColorBrewer')){
  install.packages("RColorBrewer", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}
if(!require('doParallel')){
  install.packages("doParallel", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}
if(!require('parallel')){
  install.packages("parallel", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}
if(!require('foreach')){
  install.packages("foreach", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}
if(!require('here')){
  install.packages("here", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}
if(!require('fs')){
  install.packages("fs", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}
if(!require('ggcorrplot')){
  install.packages("ggcorrplot", lib = '~/R/x86_64-pc-linux-gnu-library/3.0', repos = 'https://cran.cnr.berkeley.edu/')
}

# set lib path #
.libPaths('~/R/x86_64-pc-linux-gnu-library/3.0')

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
