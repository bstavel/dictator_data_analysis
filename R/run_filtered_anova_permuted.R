run_filtered_anova_permuted <- function(results, brain_behave_data, region_name, all_results = FALSE) {
  
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
  
  anova_chisq_dis <- NULL
  anova_chisq_adv <- NULL
  
  
  # disadvantageous anovas #
  
  if(nrow(filtered_disadvantageous) > 0) {
    for(row in 1:nrow(filtered_disadvantageous)) {
      
      bin <- as.character(filtered_disadvantageous[row, "bin"])
      elec <- as.character(filtered_disadvantageous[row, "electrode"])
      eval(parse(text = paste0("ineq_model <- lm(", bin, "~ ineq_disadvent + self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
      eval(parse(text = paste0("base_model <- lm(", bin, "~ self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
      
      anova_sum <- lrtest(base_model, ineq_model)
      anova_dis_pval[row] <- anova_sum$`Pr(>Chisq)`[2]
      anova_chisq_dis[row] <- anova_sum$Chisq[2]
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
      
      anova_sum <- lrtest(base_model, ineq_model)
      anova_adv_pval[row] <- anova_sum$`Pr(>Chisq)`[2]
      anova_chisq_adv[row] <- anova_sum$Chisq[2]
      diff_adj_r2_adv[row] <- summary(ineq_model)$adj.r.squared - summary(base_model)$adj.r.squared
      perct_adj_r2_adv[row] <- (summary(ineq_model)$adj.r.squared - summary(base_model)$adj.r.squared)/summary(base_model)$adj.r.squared
      
    }
  }
  
  sum_chisq_bins <- function(indices, values) {
    sum_temp <- values[indices[1]]
    back_idx <- indices[1]
    sum_vec <- rep(0, length(values))
    for(idx in 1:length(indices)){
      
      if(is.na(indices[idx + 1])) {
         if((indices[idx] - indices[idx-1] == 1) ) {
            sum_vec[back_idx:indices[idx]] <- sum_temp
          } else {
            sum_vec[indices[idx]] <- values[indices[idx]]
          }
        
      } else if(indices[idx + 1] - indices[idx] == 1) {
        sum_temp <- sum_temp + values[indices[idx + 1]]
        
      } else if(idx != 1){
        if(indices[idx] - indices[idx-1] == 1) {
        sum_vec[back_idx:indices[idx]] <- sum_temp
        back_idx <- indices[idx + 1]
        sum_temp <- values[indices[idx + 1]]
        } else {
          sum_vec[indices[idx]] <- values[indices[idx]]
        }
      } else {
        sum_vec[indices[idx]] <- values[indices[idx]]
      }
    }

    return(sum_vec)
  }
  
  ## Find out the longest stretch in which pval < 0.05 and create sum-of-F-stats statistic ##
  # disadvantageous #
  dis_chisq_stretch <- NULL
  indices <- which(anova_dis_pval < 0.05)
  
  if(is.na(indices[1])) {
    # if no stretch just take the max vals #
    dis_chisq_stretch[indices] <- max(anova_chisq_dis)
    
  } else {
    # if a stretch take the sum of f stats #  
    dis_chisq_stretch <-  sum_chisq_bins(indices, anova_chisq_dis)
  }
  
  # advantageous #
  adv_chisq_stretch <- NULL
  indices <- which(anova_adv_pval < 0.05)
  if(is.na(indices[1])) {
    # if no stretch just take the max vals #
    adv_chisq_stretch[indices] <- max(anova_chisq_adv)
    
  } else {
    # if a stretch take the sum of f stats #  
    adv_chisq_stretch <-  sum_chisq_bins(indices, anova_chisq_adv) # Summary stat
  }
    
  filtered_disadvantageous$anova_p <- anova_dis_pval
  filtered_advantageous$anova_p <- anova_adv_pval
  filtered_disadvantageous$diff_adj_r2 <- diff_adj_r2_dis
  filtered_advantageous$diff_adj_r2 <- diff_adj_r2_adv
  filtered_disadvantageous$perct_adj_r2 <- perct_adj_r2_dis
  filtered_advantageous$perct_adj_r2 <- perct_adj_r2_adv
  filtered_disadvantageous$chisq_stretch <- dis_chisq_stretch
  filtered_advantageous$chisq_stretch <- adv_chisq_stretch
  
  filtered_advantageous <- filtered_advantageous %>%   mutate(anova_permuted = 0)
  for(elec in unique(results$electrode)){  
    filtered_advantageous <- filtered_advantageous %>%
      group_by(chisq_stretch) %>%
      mutate_cond(electrode == elec, anova_permuted = (sum(chisq_stretch < null_presentation_data_adv[elec, ]))/1000 )
    
  } 
  filtered_disadvantageous <- filtered_disadvantageous %>%   mutate(anova_permuted = 0)
  for(elec in unique(results$electrode)){  
    filtered_disadvantageous <- filtered_disadvantageous %>%
      group_by(chisq_stretch) %>%
      mutate_cond(electrode == elec, anova_permuted = (sum(chisq_stretch < null_presentation_data_dis[elec, ]))/1000 )
    
  } 
  # save results to results folder #
  write.csv(filtered_disadvantageous, path(here(), "results", paste0(region_name, "_", tag, "_anova_results_disadvantageous_permuted.csv")))
  write.csv(filtered_advantageous, path(here(), "results", paste0(region_name, "_", tag, "_anova_results_advantageous_permuted.csv")))
  
}
