run_permuted_multiple_regressions_par <- function(brain_behave_data, electrodes, regressor, nBins, region_name, niter = 1000, sub, tag) {
  print(paste0("Beginning regressions for predictor: ", regressor))
  for(elec in electrodes) {
    
    # update #
    print(paste0("Starting electrode ", which(electrodes %in% elec), " of ", length(electrodes)))
    
    # create results dfs #
    results <- data.frame(matrix(nrow = length(nBins), ncol = 10))
    colnames(results) <- c("electrode", "bin", "R2", "Beta_1", "Beta_2", "Fstat", "p_1", "p_2", "fstretch", "perm_p")
    
    # fill results dfs with electrode & bin info #
    results$electrode <- elec
    results$bin <- nBins
    
    # filter to single electrode #
    brain_behave_data_elec <- brain_behave_data %>% filter(electrodes == elec)
    # brain_behave_data_elec <- brain_behave_data_elec[complete.cases(brain_behave_data_elec), ]
    
    # initialize temp vars #
    reg_1 <- regressor[1]
    reg_2 <- regressor[2]
    r2 <- NULL
    fstat <- NULL
    beta_1 <- NULL
    beta_2 <- NULL
    lm_pval_1 <- NULL
    lm_pval_2 <- NULL
    
    for (bin in nBins) {
      # run models #
      bin_vec <- brain_behave_data_elec[, bin]
      reg_1_vec <- brain_behave_data_elec[, reg_1]
      reg_2_vec <- brain_behave_data_elec[, reg_2]
      model <- summary(lm(bin_vec ~ reg_1_vec + reg_2_vec))
      
      # store info from models #
      r2[bin] <- model$r.squared
      fstat[bin] <- model$fstatistic[1]
      beta_1[bin] <- model$coefficients[2, 1]
      beta_2[bin] <- model$coefficients[3, 1]
      lm_pval_1[bin] <- model$coefficients[2,4]
      lm_pval_2[bin] <- model$coefficients[3,4]
      
    }
    
    ## Find out the longest stretch in which pval < 0.05 and create sum-of-F-stats statistic ##
    stretch <- lm_pval_1 < 0.05
    indices <- stretch_start_end(stretch)
    
    # if there is no stretch, take the max beta and max f statistic #
    if(is.na(indices[1])) {
      beta_stretch <- max(beta_1)
      fstat_stretch <- max(fstat)
      # if there is a stretch, sum the betas and the f statistics #
    } else {
      beta_stretch <- sum(beta_1[indices[1]:indices[2]]) # Summary stat
      fstat_stretch <- sum(fstat[indices[1]:indices[2]]) # Summary stat
    }
    
    # only run permutation testing if fstratch is above 7, otherwise extremeley unlikely to be sig
    if( fstat_stretch > 1 ) {
      # Create null distribution by shuffling labels #
      null_fstat_stretch <- NULL
      null_fstat_stretch <- foreach(h = 1:niter, .inorder=FALSE) %dopar% {
        # print status #
        if(h %% 100 == 0){
          print(paste0( (h/niter) * 100, "% complete for electrode", which(electrodes %in% elec), " of ", length(electrodes)))
        }
        
        # intialize vars #
        null_fstat <- NULL
        null_lm_pval <- NULL
        
        for (bin in nBins) { # This is the slooooow step
          
          # run models #
          bin_vec <- brain_behave_data_elec[, bin]
          reg_1_vec <- brain_behave_data_elec[, reg_1]
          reg_2_vec <- brain_behave_data_elec[, reg_2]
          null_lm <-summary(lm(bin_vec ~ sample(reg_1_vec) + reg_1_vec, data = brain_behave_data_elec))
          
          # save info from models #
          null_fstat[bin] <- null_lm$fstatistic[1]
          null_lm_pval[bin] <- null_lm$coefficients[2,4]
          
        }
        
        ## Find out the longest stretch in which pval < 0.05 and create sum-of-F-stats statistic ##
        # disadvantageous #
        stretch <- null_lm_pval < 0.05
        indices <- stretch_start_end(stretch)
        if(is.na(indices[1])) {
          # if no stretch just take the max vals #
          null_fstat_stretch[h] <- max(null_fstat)
          
        } else {
          # if a stretch take the sum of f stats #  
          null_fstat_stretch[h] <- sum(null_fstat[indices[1]:indices[2]]) # Summary stat
        }
        
        return(null_fstat_stretch[h])
      }
      
      # our test statistic for the permutation test is the sum of f tests #
      perm_pval <- sum(null_fstat_stretch > fstat_stretch)/niter
      
    } else {
      perm_pval <- 2
    }
    
    # save vals in results df #
    results <- results %>% 
      mutate_cond(electrode == elec, R2 = r2, Beta_1 = beta_1, Beta_2 = beta_2, Fstat = fstat, p_1 = lm_pval_1, p_2 = lm_pval_2, fstretch = fstat_stretch, perm_p = perm_pval)
    
    
    # save results to results folder #
    write.csv(results, path(here(), "results", sub, paste0(region_name, "_", elec, "_", reg_1, "_", reg_2, "_", tag, "_results.csv")))
    
  }
  
}
