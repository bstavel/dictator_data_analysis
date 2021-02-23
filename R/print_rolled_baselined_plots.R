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
library(viridis)
library(beepr)
library(lmtest)

## hand written functions ##
source(path(here(), "R", "load_behave_data.R"))
source(path(here(), "R", "prep_behave_data.R"))
source(path(here(), "R", "stretch_start_end.R"))
source(path(here(), "R", "load_high_gamma_data.R"))
source(path(here(), "R", "rolling_window_and_baseline.R"))
source(path(here(), "R", "run_permuted_regressions.R"))
source(path(here(), "R", "compile_results.R"))
source(path(here(), "R", "run_filtered_anova.R"))
source(path(here(), "R", "merge_stepped_anova_results.R"))
source(path(here(), "R", 'mutate_cond.R'))

## plotting helpers ##
# ggthemr("solarized")

# ## paralellization ##
# nCores <- 2
# registerDoParallel(nCores)

## function to filter to ROI, baseline, compress, and c=merge within behvaioral data ##
print_IR35_presentation_visualization_data <- function(file_path_to_power_data, clean_name, band, sub) {
  
  # necessary paths #
  file_path_to_electrode_names <- path(here(), "munge", "savio_cluster", paste0(sub, "_electrodes_presentation_locked_rscaler_2575.csv"))
  file_path_to_elecs_of_interest <- path(here(), "munge", paste0(sub, "_elecs_of_interest.csv"))
  
  # load data #
  hg_data <- load_high_gamma_data(file_path_to_power_data, file_path_to_electrode_names)
  elecs_to_use <- read.csv(file_path_to_elecs_of_interest)
  
  # filter to good electrodes #
  hg_data <- hg_data %>% filter(grepl(paste(elecs_to_use$Electrode, collapse = "|"), electrodes))
  
  # average and baseline #
  hg_clean_rolled <- rolling_window_and_baseline(hg_data, lWin = 200, lOver = 50, choice_locked = F)
  
  plot_with_base_and_roll <- hg_clean_rolled %>%
    gather(key = "bin", value = "power", -electrodes, -trial) %>%
    group_by(electrodes, bin) %>%
    mutate(mean_power = mean(power))  %>%
    mutate(lower_sem_power = mean(power) - 2*(sd(power)/n())) %>%
    mutate(upper_sem_power = mean(power) + 2*(sd(power)/n())) %>%
    mutate(time = (50 * as.numeric(gsub("bin_", "", bin))) - 200) %>%
    ggplot(., aes(x = time, y = mean_power)) +
    geom_line() +
    geom_ribbon(aes(ymin = lower_sem_power, ymax = upper_sem_power), alpha = .5) +
    geom_vline(xintercept = 0, color = "black") +
    theme(panel.background = element_rect(fill = "white")) +
    facet_wrap(~electrodes)
  ggtitle(paste0(band, " ~ With Baseline, with Rolling, "))
  
  ggsave(path(here(), "figures", "data_checks", paste0(sub, "_", band, "_rolled_power_envelop.pdf")), plot = plot_with_base_and_roll, height = 18, width = 42, units = "in")
  
}


## Extract for all ##
subs <- c("IR26", "IR28", "IR35", "IR57")
for( sub in subs) {
  
  # delta #
  band <- "delta"
  file_path_to_power_data <- path(here(), "munge", "savio_cluster", paste0(sub, "_delta_munge_presentation_locked_rscaler_2575.csv"))
  print_IR35_presentation_visualization_data(file_path_to_power_data, clean_name, band, sub)
  # theta #
  band <- "theta"
  file_path_to_power_data <- path(here(), "munge", "savio_cluster", paste0(sub, "_theta_munge_presentation_locked_rscaler_2575.csv"))
  print_IR35_presentation_visualization_data(file_path_to_power_data, clean_name, band, sub)
  # alpha #
  band <- "alpha"
  file_path_to_power_data <- path(here(), "munge", "savio_cluster", paste0(sub, "_alpha_munge_presentation_locked_rscaler_2575.csv"))
  print_IR35_presentation_visualization_data(file_path_to_power_data, clean_name, band, sub)
  # beta #
  band <- "beta"
  file_path_to_power_data <- path(here(), "munge", "savio_cluster", paste0(sub, "_beta_munge_presentation_locked_rscaler_2575.csv"))
  print_IR35_presentation_visualization_data(file_path_to_power_data, clean_name, band, sub)
  # gamma #
  band <- "gamma"
  file_path_to_power_data <- path(here(), "munge", "savio_cluster", paste0(sub, "_gamma_munge_presentation_locked_rscaler_2575.csv"))
  print_IR35_presentation_visualization_data(file_path_to_power_data, clean_name, band, sub)
  # hfa #
  band <- "hfa"
  file_path_to_power_data <- path(here(), "munge", "savio_cluster", paste0(sub, "_hfa_munge_presentation_locked_rscaler_2575.csv"))
  print_IR35_presentation_visualization_data(file_path_to_power_data, clean_name, band, sub)
  
  
}
