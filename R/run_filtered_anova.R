run_filtered_anova <- function(results, brain_behave_data, sub, region_name, all_results = FALSE, type) {
  
  if(all_results == F){
    # get active electrodes #
    filtered_disadvantageous <- results %>% filter(predictor == "ineq_disadvent") %>% filter(perm_p < 0.05 & p < 0.05)
    filtered_advantageous <-  results %>% filter(predictor == "ineq_advent") %>% filter(perm_p < 0.05 & p < 0.05)
    dis_correction <- nrow(filtered_disadvantageous)
    adv_correction <- nrow(filtered_advantageous)
  } else  {
    # get all electrodes #
    filtered_disadvantageous <- results %>% filter(predictor == "ineq_disadvent")
    filtered_advantageous <-  results %>% filter(predictor == "ineq_advent") 
    dis_correction <- nrow(filtered_disadvantageous %>% filter(perm_p < 0.05 & p < 0.05))
    adv_correction <- nrow(filtered_advantageous %>% filter(perm_p < 0.05 & p < 0.05))
  }
  
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
      eval(parse(text = paste0("ineq_model <- lm(", bin, "~ ineq_advent + self_var_payoff + other_var_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
      eval(parse(text = paste0("base_model <- lm(", bin, "~ self_var_payoff + other_var_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
      
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
  filtered_disadvantageous$correction <- dis_correction
  filtered_advantageous$correction <- adv_correction
  
  # save results to results folder #
  write.csv(filtered_disadvantageous, path(here(), "results", sub, paste0(region_name, "_", type, "_anova_results_disadvantageous_test2.csv")))
  write.csv(filtered_advantageous, path(here(), "results", sub, paste0(region_name, "_", type, "_anova_results_advantageous_test2.csv")))
  
}
