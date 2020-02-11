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
      
      # Create null distribution by shuffling labels #
      null_chisq_stretch <- NULL
      null_chisq <- NULL
      null_anova_pval <- NULL
      null_chisq_stretch <- foreach(h = 1:niter, .inorder=FALSE) %dopar% {
        # print status #
        if(h %% 100 == 0){
          print(paste0( (h/niter) * 100, "% complete for electrode", which(electrodes %in% elec), " of ", length(electrodes)))
        }
        

        # run models #
        set.seed(h)
        eval(parse(text = paste0("ineq_model <- lm(", bin, "~ sample(ineq_disadvent) + self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
        eval(parse(text = paste0("base_model <- lm(", bin, "~ self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
        null_anova_sum <- lrtest(base_model, ineq_model)

        
        # save info from models #
        null_chisq[row] <- null_anova_sum$Chisq[2]
        null_anova_pval[row] <- null_anova_sum$`Pr(>Chisq)`[2]

        
        ## Find out the longest stretch in which pval < 0.05 and create sum-of-F-stats statistic ##
        # disadvantageous #
        stretch <- null_anova_pval < 0.05
        indices <- stretch_start_end(stretch)
        if(is.na(indices[1])) {
          # if no stretch just take the max vals #
          null_chisq_stretch[h] <- max(null_chisq)
          
        } else {
          # if a stretch take the sum of f stats #  
          null_chisq_stretch[h] <- sum(null_chisq[indices[1]:indices[2]]) # Summary stat
        }
        
        return(null_chisq_stretch[h])
      }
      
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
      anova_adv_pval[row] <- anova_sum$`Pr(>Chisq)`[2]
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
  write.csv(filtered_disadvantageous, path(here(), "results", paste0(region_name, "_anova_results_disadvantageous.csv")))
  write.csv(filtered_advantageous, path(here(), "results", paste0(region_name, "_anova_results_advantageous.csv")))
  
}
