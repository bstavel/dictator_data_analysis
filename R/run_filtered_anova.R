run_filtered_anova <- function(results, brain_behave_data, region_name) {
  
  # get active electrodes #
  filtered_disadvantageous <- results %>% filter(predictor == "ineq_disadvent") %>% filter(perm_p < 0.05 & p < 0.05)
  filtered_advantageous <-  results %>% filter(predictor == "ineq_advent") %>% filter(perm_p < 0.05 & p < 0.05)
  
  # initialize values #
  anova_dis_pval <- NULL
  anova_adv_pval <- NULL
  
  diff_adj_r2_dis <- NULL
  diff_adj_r2_adv <- NULL

  perct_adj_r2_dis <- NULL  
  perct_adj_r2_adv <- NULL  
  
  # disadvantageous anovas #
  
  if(nrow(filtered_disadvantageous) > 0) {
    for(row in 1:nrow(filtered_disadvantageous)) {
      
      bin <- as.character(filtered_disadvantageous[row, "bin"])
      elec <- as.character(filtered_disadvantageous[row, "electrode"])
      eval(parse(text = paste0("ineq_model <- lm(", bin, "~ ineq_disadvent + self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
      eval(parse(text = paste0("base_model <- lm(", bin, "~ self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
      
      anova_sum <- anova(base_model, ineq_model, test = "Chisq")
      anova_dis_pval[row] <- anova_sum$`Pr(>Chi)`[2]
      diff_adj_r2_dis[row] <- summary(ineq_model)$adj.r.squared - summary(base_model)$adj.r.squared
      perct_adj_r2_dis[row] <- (summary(ineq_model)$adj.r.squared - summary(base_model)$adj.r.squared)/summary(base_model)$adj.r.squared
      
    }
  }
  
  # advantageous anovas #
  if(nrow(filtered_advantageous) > 0) {
    for(row in 1:nrow(filtered_advantageous)) {
      
      bin <- as.character(filtered_advantageous[row, "bin"])
      elec <- as.character(filtered_advantageous[row, "electrode"])
      eval(parse(text = paste0("ineq_model <- lm(", bin, "~ ineq_advent + self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
      eval(parse(text = paste0("base_model <- lm(", bin, "~ self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
      
      anova_sum <- anova(base_model, ineq_model, test = "Chisq")
      anova_adv_pval[row] <- anova_sum$`Pr(>Chi)`[2]
      diff_adj_r2_adv[row] <- summary(ineq_model)$adj.r.squared - summary(base_model)$adj.r.squared
      perct_adj_r2_adv[row] <- (summary(ineq_model)$adj.r.squared - summary(base_model)$adj.r.squared)/summary(base_model)$adj.r.squared
      
    }
  }
    
  filtered_disadvantageous$anova_p <- anova_dis_pval
  filtered_advantageous$anova_p <- anova_adv_pval
  filtered_disadvantageous$diff_adj_r2 <- diff_adj_r2_dis
  filtered_advantageous$diff_adj_r2 <- diff_adj_r2_adv
  filtered_disadvantageous$perct_adj_r2 <- perct_adj_r2_dis
  filtered_advantageous$perct_adj_r2 <- perct_adj_r2_adv
  
  # save results to results folder #
  write.csv(filtered_disadvantageous, path(here(), "results", paste0(region_name, "_anova_results_disadvantageous.csv")))
  write.csv(filtered_advantageous, path(here(), "results", paste0(region_name, "_anova_results_advantageous.csv")))
  
}
