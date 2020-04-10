run_permuted_regressions <- function(brain_behave_data, electrodes, regressor, nBins, region_name, niter = 1000) {
  
  foreach(elecIdx = 1:length(electrodes))  %dopar% {
    # get elec #
    elec <- electrodes[elecIdx]
    
    # create results dfs #
    results <- data.frame(matrix(nrow = length(nBins), ncol = 8))
    colnames(results) <- c("electrode", "bin", "R2", "beta", "Fstat", "p", "fstretch", "perm_p")
    
    # fill results dfs with electrode & bin info #
    results$electrode <- elec
    results$bin <- paste0("bin_", 1:nBins)
    
    # filter to single electrode #
    brain_behave_data_elec <- brain_behave_data %>% filter(electrodes == elec)
    # brain_behave_data_elec <- brain_behave_data_elec[complete.cases(brain_behave_data_elec), ]
    
    # initialize temp vars #
    r2 <- NULL
    fstat <- NULL
    beta <- NULL
    lm_pval <- NULL
    
    for (bin in nBins) {
      # run models #
      bin_vec <- brain_behave_data_elec[, bin]
      reg_vec <- brain_behave_data_elec[, regressor]
      model <- summary(lm(bin_vec ~ reg_vec))
      
      # store info from models #
      r2[bin] <- model$r.squared
      fstat[bin] <- model$fstatistic[1]
      beta[bin] <- model$coefficients[2, 1]
      lm_pval[bin] <- model$coefficients[2,4]
      
    }
    
    ## Find out the longest stretch in which pval < 0.05 and create sum-of-F-stats statistic ##
    stretch <- lm_pval < 0.05
    indices <- stretch_start_end(stretch)  # what does stretch_start_end do???
    
    # if there is no stretch, take the max beta and max f statistic #
    if(is.na(indices[1])) {
      fstat_stretch <- max(fstat)
      # if there is a stretch, sum the betas and the f statistics #
    } else {
      fstat_stretch <- sum(fstat[indices[1]:indices[2]]) # Summary stat
    }
    
    
    # Create null distribution by shuffling labels #
    null_fstat_stretch <- NULL
    for (h in 1:niter) {
      # intialize vars #
      null_fstat <- NULL
      null_lm_pval <- NULL
      
     for (bin in nBins) { # This is the slooooow step
        
        # if(remove_outliers) {
        #   temp <- data.frame(peak=peak)
        #   temp <- remove_outliers(temp,peak)
        #   peak <- temp$peak
        # }
        
        # run models #
        bin_vec <- brain_behave_data_elec[, bin]
        reg_vec <- brain_behave_data_elec[, regressor]
        null_lm <-summary(lm(bin_vec ~ sample(reg_vec), data = brain_behave_data_elec))
        
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
      
    }
    
    # our test statistic for the permutation test is the sum of f tests #
    perm_pval <- sum(null_fstat_stretch > fstat_stretch)/niter

    # save vals in results df #
    results <- results %>% 
      mutate_cond(electrode == elec, R2 = r2, beta = abs_beta, Fstat = fstat, p = lm_pval, fstretch = fstat_stretch, perm_p = perm_pval)
    
  
    # save results to results folder #
    write.csv(results, path(here(), "results", paste0(region_name, "_", elec, "_", regressor, "_results.csv")))
  
  }
  
}
