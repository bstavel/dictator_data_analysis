

## get null distributions ##
niter <- 1000
nCores <- 2
registerDoParallel(nCores)

# presesntation #
path_hg_clean <- path(here(), "munge", "hg_behave_choice_locked_cut_fixation.csv")
hg_behave <-  read.csv(path_hg_clean)
results <- choice_results
brain_behave_data <- hg_behave
region_name <- "All"
all_results = TRUE

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

null_chisq_stretch_dis <- matrix(nrow= length(unique(results$electrode)), ncol = niter)
row.names(null_chisq_stretch_dis) <- unique(results$electrode)
# adv
null_chisq_stretch_adv <- matrix(nrow= length(unique(results$electrode)), ncol = niter)
row.names(null_chisq_stretch_adv) <- unique(results$electrode)
for(elec in unique(results$electrode)) {
  null_chisq_stretch_temp <- NULL
  null_chisq_stretch_temp_adv <- NULL
  null_chisq_stretch_temp_dis <- NULL
  
  # run par loop over 1000 samples #
  null_chisq_stretch_temp <- foreach(h = 1:niter, .inorder=FALSE, .combine = cbind) %dopar% {
    # print status #
    if(h %% 100 == 0){
      print(paste0( (h/niter) * 100, "% complete"))
    }
    
    null_chisq_adv <- NULL
    null_anova_pval_adv <- NULL      
    null_chisq_dis <- NULL
    null_anova_pval_dis <- NULL
    for(bin in unique(results$bin)){
      # run models #
      set.seed(h)
      eval(parse(text = paste0("ineq_model <- lm(", bin, "~ sample(ineq_advent) + self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
      eval(parse(text = paste0("base_model <- lm(", bin, "~ self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
      null_anova_sum <- lrtest(base_model, ineq_model)
      
      # save info from models #
      null_chisq_adv[bin] <- null_anova_sum$Chisq[2]
      null_anova_pval_adv[bin] <- null_anova_sum$`Pr(>Chisq)`[2]
      
      set.seed(h)
      eval(parse(text = paste0("ineq_model <- lm(", bin, "~ sample(ineq_disadvent) + self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
      eval(parse(text = paste0("base_model <- lm(", bin, "~ self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
      null_anova_sum <- lrtest(base_model, ineq_model)        
      
      # save info from models #
      null_chisq_dis[bin] <- null_anova_sum$Chisq[2]
      null_anova_pval_dis[bin] <- null_anova_sum$`Pr(>Chisq)`[2]
      
    }
    
    ## Find out the longest stretch in which pval < 0.05 and create sum-of-F-stats statistic ##
    # disadvantageous #
    stretch <- null_anova_pval_dis < 0.05
    indices <- stretch_start_end(stretch)
    if(is.na(indices[1])) {
      # if no stretch just take the max vals #
      null_chisq_stretch_temp_dis[h] <- max(null_chisq_dis)
      
    } else {
      # if a stretch take the sum of f stats #  
      null_chisq_stretch_temp_dis[h] <- sum(null_chisq_dis[indices[1]:indices[2]]) # Summary stat
    }
    
    # advantageous #
    stretch <- null_anova_pval_adv < 0.05
    indices <- stretch_start_end(stretch)
    if(is.na(indices[1]) | indices[1] == indices[2]) {
      # if no stretch just take the max vals #
      null_chisq_stretch_temp_adv[h] <- max(null_chisq_adv)
      
    } else {
      # if a stretch take the sum of f stats #  
      null_chisq_stretch_temp_adv[h] <- sum(null_chisq_adv[indices[1]:indices[2]]) # Summary stat
    }      
    
    null_chisq_stretch_temp <- rbind(null_chisq_stretch_temp_adv[h], null_chisq_stretch_temp_dis[h])
    return(null_chisq_stretch_temp)
  }
  
  # update final matrices #
  null_chisq_stretch_adv[elec, 1:niter] <- null_chisq_stretch_temp[1, ]
  null_chisq_stretch_dis[elec, 1:niter] <- null_chisq_stretch_temp[2, ]
  
}

null_presentation_data_adv <- null_chisq_stretch_adv
null_presentation_data_dis <- null_chisq_stretch_dis


## choice ##
path_hg_clean <- path(here(), "munge", "hg_behave_presentation_locked.csv")
hg_behave <-  read.csv(path_hg_clean)
results <- presentation_results
brain_behave_data <- hg_behave
region_name <- "All"
all_results = TRUE

null_chisq_stretch_dis <- matrix(nrow= length(unique(results$electrode)), ncol = niter)
row.names(null_chisq_stretch_dis) <- unique(results$electrode)
# adv
null_chisq_stretch_adv <- matrix(nrow= length(unique(results$electrode)), ncol = niter)
row.names(null_chisq_stretch_adv) <- unique(results$electrode)
for(elec in unique(results$electrode)) {
  null_chisq_stretch_temp <- NULL
  null_chisq_stretch_temp_adv <- NULL
  null_chisq_stretch_temp_dis <- NULL
  
  # run par loop over 1000 samples #
  null_chisq_stretch_temp <- foreach(h = 1:niter, .inorder=FALSE, .combine = cbind) %dopar% {
    # print status #
    if(h %% 100 == 0){
      print(paste0( (h/niter) * 100, "% complete"))
    }
    
    null_chisq_adv <- NULL
    null_anova_pval_adv <- NULL      
    null_chisq_dis <- NULL
    null_anova_pval_dis <- NULL
    for(bin in unique(results$bin)){
      # run models #
      set.seed(h)
      eval(parse(text = paste0("ineq_model <- lm(", bin, "~ sample(ineq_advent) + self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
      eval(parse(text = paste0("base_model <- lm(", bin, "~ self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
      null_anova_sum <- lrtest(base_model, ineq_model)
      
      # save info from models #
      null_chisq_adv[bin] <- null_anova_sum$Chisq[2]
      null_anova_pval_adv[bin] <- null_anova_sum$`Pr(>Chisq)`[2]
      
      set.seed(h)
      eval(parse(text = paste0("ineq_model <- lm(", bin, "~ sample(ineq_disadvent) + self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
      eval(parse(text = paste0("base_model <- lm(", bin, "~ self_payoff + other_payoff, data = brain_behave_data[brain_behave_data$electrodes == '", elec, "', ])")))
      null_anova_sum <- lrtest(base_model, ineq_model)        
      
      # save info from models #
      null_chisq_dis[bin] <- null_anova_sum$Chisq[2]
      null_anova_pval_dis[bin] <- null_anova_sum$`Pr(>Chisq)`[2]
      
    }
    
    ## Find out the longest stretch in which pval < 0.05 and create sum-of-F-stats statistic ##
    # disadvantageous #
    stretch <- null_anova_pval_dis < 0.05
    indices <- stretch_start_end(stretch)
    if(is.na(indices[1]) | indices[1] == indices[2]) {
      # if no stretch just take the max vals #
      null_chisq_stretch_temp_dis[h] <- max(null_chisq_dis)
      
    } else {
      # if a stretch take the sum of f stats #  
      null_chisq_stretch_temp_dis[h] <- sum(null_chisq_dis[indices[1]:indices[2]]) # Summary stat
    }
    
    # advantageous #
    stretch <- null_anova_pval_adv < 0.05
    indices <- stretch_start_end(stretch)
    if(is.na(indices[1]) | indices[1] == indices[2]) {
      # if no stretch just take the max vals #
      null_chisq_stretch_temp_adv[h] <- max(null_chisq_adv)
      
    } else {
      # if a stretch take the sum of f stats #  
      null_chisq_stretch_temp_adv[h] <- sum(null_chisq_adv[indices[1]:indices[2]]) # Summary stat
    }      
    
    null_chisq_stretch_temp <- rbind(null_chisq_stretch_temp_adv[h], null_chisq_stretch_temp_dis[h])
    return(null_chisq_stretch_temp)
  }
  
  # update final matrices #
  null_chisq_stretch_adv[elec, 1:niter] <- null_chisq_stretch_temp[1, ]
  null_chisq_stretch_dis[elec, 1:niter] <- null_chisq_stretch_temp[2, ]
  
}


null_choice_data_adv <- null_chisq_stretch_adv
null_choice_data_dis <- null_chisq_stretch_dis
