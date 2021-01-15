run_permuted_interaction_regressions_par <- function(brain_behave_data, electrodes, regressor, nBins, region_name, niter = 1000, sub, tag) {
  print(paste0("Beginning regressions for predictor: ", regressor))
  for(elec in electrodes) {
    
    # update #
    print(paste0("Starting electrode ", which(electrodes %in% elec), " of ", length(electrodes)))
    
    # create results dfs #
    results <- data.frame(matrix(nrow = length(nBins), ncol = 11))
    colnames(results) <- c("electrode", "bin", "Beta_1", "Beta_2", "Beta_Inter", "Chi_Sq", "p_1", "p_2", "p_inter", "cstretch", "perm_p")
    
    # fill results dfs with electrode & bin info #
    results$electrode <- elec
    results$bin <- nBins
    
    # filter to single electrode #
    brain_behave_data_elec <- brain_behave_data %>% filter(electrodes == elec)
    # brain_behave_data_elec <- brain_behave_data_elec[complete.cases(brain_behave_data_elec), ]
    
    # initialize temp vars #
    reg_1 <- regressor[1]
    reg_2 <- regressor[2]
    chisq <- NULL
    beta_1 <- NULL
    beta_2 <- NULL
    beta_inter <- NULL
    lm_pval_1 <- NULL
    lm_pval_2 <- NULL
    lm_pval_inter <- NULL
    model_com_pval <- NULL
    
    for (bin in nBins) {
      # run models #
      bin_vec <- brain_behave_data_elec[, bin]
      reg_1_vec <- brain_behave_data_elec[, reg_1]
      reg_2_vec <- brain_behave_data_elec[, reg_2]
      base_model <- lm(bin_vec ~ reg_1_vec + reg_2_vec)
      interaction_model <- lm(bin_vec ~ reg_1_vec + reg_2_vec + reg_1_vec*reg_2_vec)
      model_comp <- lrtest(base_model, interaction_model)
      
      # store info from models #
      chisq[bin] <- model_comp$Chisq[2]
      beta_1[bin] <- summary(interaction_model)$coefficients[2, 1]
      beta_2[bin] <- summary(interaction_model)$coefficients[3, 1]
      beta_inter[bin] <-  summary(interaction_model)$coefficients[3, 1]
      lm_pval_1[bin] <- summary(interaction_model)$coefficients[2,4]
      lm_pval_2[bin] <- summary(interaction_model)$coefficients[3,4]
      lm_pval_inter[bin] <- summary(interaction_model)$coefficients[3,4]
      model_com_pval[bin] <-  model_comp$`Pr(>Chisq)`[2]
      
    }
    
    ## Find out the longest stretch in which pval < 0.05 and create sum-of-F-stats statistic ##
    stretch <- model_com_pval < 0.05
    indices <- stretch_start_end(stretch)
    
    # if there is no stretch, take the max beta and max f statistic #
    if(is.na(indices[1])) {
      chisq_stretch <- max(chisq)
      # if there is a stretch, sum the betas and the f statistics #
    } else {
      chisq_stretch <- sum(chisq[indices[1]:indices[2]]) # Summary stat
    }
    
    # only run permutation testing if fstratch is above 7, otherwise extremeley unlikely to be sig
    if( chisq_stretch > 5 ) {
      # Create null distribution by shuffling labels #
      null_chisq_stretch <- NULL
      null_chisq_stretch <- foreach(h = 1:niter, .inorder=FALSE) %dopar% {
        # print status #
        if(h %% 100 == 0){
          print(paste0( (h/niter) * 100, "% complete for electrode", which(electrodes %in% elec), " of ", length(electrodes)))
        }
        
        # intialize vars #
        null_chisq <- NULL
        null_model_com_pval <- NULL
        
        for (bin in nBins) { # This is the slooooow step
          
          # run models #
          bin_vec <- brain_behave_data_elec[, bin]
          reg_1_vec <- brain_behave_data_elec[, reg_1]
          reg_2_vec <- brain_behave_data_elec[, reg_2]
          randomization <- sample(1:length(reg_1_vec))
          
          null_base_model <- lm(bin_vec ~ reg_1_vec[randomization] + reg_2_vec[randomization])
          null_interaction_model <- lm(bin_vec ~ reg_1_vec[randomization] + reg_2_vec[randomization] + reg_1_vec[randomization]*reg_2_vec[randomization])
          null_model_comp <- lrtest(base_model, interaction_model)
          
          # save info from models #
          null_chisq[bin] <- null_model_comp$Chisq[2]
          null_model_com_pval[bin] <- null_model_comp$`Pr(>Chisq)`[2]
          
        }
        
        ## Find out the longest stretch in which pval < 0.05 and create sum-of-F-stats statistic ##
        # disadvantageous #
        stretch <- null_model_com_pval < 0.05
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
      
      # our test statistic for the permutation test is the sum of f tests #
      perm_pval <- sum(null_chisq_stretch > chisq_stretch)/niter
      
    } else {
      perm_pval <- 2
    }
    
    # save vals in results df #
    results <- results %>% 
      mutate(Beta_1 = beta_1, Beta_2 = beta_2, Beta_Inter = beta_inter, Chi_Sq = chisq, p_1 = lm_pval_1, p_2 = lm_pval_2, p_inter = lm_pval_inter, cstretch = chisq_stretch, perm_p = perm_pval)
    
    
    # save results to results folder #
    write.csv(results, path(here(), "results", sub, paste0(region_name, "_", elec, "_", reg_1, "_", reg_2, "_", tag, "_results.csv")))
    
  }
  
}
