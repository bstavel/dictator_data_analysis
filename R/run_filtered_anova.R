run_filtered_anova <- function(brain_behave_data, region_name) {
  
  # read back in the data #
  results_disadvantageous <- read.csv(path(here(), "results", paste0(region_name, "_results_disadvantageous.csv")))
  results_advantageous <- read.csv(path(here(), "results", paste0(region_name, "_results_advantageous.csv")))
  
  # get active electrodes #
  filtered_disadvantageous <- results_disadvantageous %>% filter(perm_p < 0.05 & p < 0.05)
  filtered_advantageous <- results_advantageous %>% filter(perm_p < 0.05 & p < 0.05)
  
  # initialize values #
  anova_dis_pval <- NULL
  anova_adv_pval <- NULL
  
  # disadvantageous anovas #
  for(row in 1:nrow(filtered_disadvantageous)) {
    
    bin <- as.character(filtered_disadvantageous[row, "bin"])
    elec <- as.character(filtered_disadvantageous[row, "electrode"])
    eval(parse(text = paste0("ineq_model <- lm(", bin, "~ ineq_disadvent + self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
    eval(parse(text = paste0("base_model <- lm(", bin, "~ self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
    
    anova_sum <- anova(base_model, ineq_model, test = "Chisq")
    anova_dis_pval[row] <- anova_sum$`Pr(>Chi)`[2]
    
  }
  
  # advantageous anovas #
  for(row in 1:nrow(filtered_advantageous)) {
    
    bin <- as.character(filtered_advantageous[row, "bin"])
    elec <- as.character(filtered_advantageous[row, "electrode"])
    eval(parse(text = paste0("ineq_model <- lm(", bin, "~ ineq_advent + self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
    eval(parse(text = paste0("base_model <- lm(", bin, "~ self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
    
    anova_sum <- anova(base_model, ineq_model, test = "Chisq")
    anova_adv_pval[row] <- anova_sum$`Pr(>Chi)`[2]
    
  }
  
  filtered_disadvantageous$anova_p <- anova_dis_pval
  filtered_advantageous$anova_p <- anova_adv_pval
  
  # save results to results folder #
  write.csv(filtered_disadvantageous, path(here(), "results", paste0(region_name, "_anova_results_disadvantageous.csv")))
  write.csv(filtered_advantageous, path(here(), "results", paste0(region_name, "_anova_results_advantageous.csv")))
  
}
