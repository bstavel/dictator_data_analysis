run_permuted_regressions_par <- function(brain_behave_data, electrodes, regressor, nBins, region_name, niter = 1000, sub, tag) {
  print(paste0("Beginning regressions for predictor: ", regressor))
  for(elec in electrodes) {

    # update #
    print(paste0("Starting electrode ", which(electrodes %in% elec), " of ", length(electrodes)))
    
    # create results dfs #
    results <- data.frame(matrix(nrow = length(nBins), ncol = 8))
    colnames(results) <- c("electrode", "bin", "R2", "absBeta", "Fstat", "p", "fstretch", "perm_p")
    
    # fill results dfs with electrode & bin info #
    results$electrode <- elec
    results$bin <- nBins
    
    # filter to single electrode #
    brain_behave_data_elec <- brain_behave_data %>% filter(electrodes == elec)
    # brain_behave_data_elec <- brain_behave_data_elec[complete.cases(brain_behave_data_elec), ]
    
    # initialize temp vars #
    r2 <- NULL
    fstat <- NULL
    abs_beta <- NULL
    lm_pval <- NULL
    
    for (bin in nBins) {
      # run models #
      bin_vec <- brain_behave_data_elec[, bin]
      reg_vec <- brain_behave_data_elec[, regressor]
      model <- summary(lm(bin_vec ~ reg_vec))
      
      # store info from models #
      r2[bin] <- model$r.squared
      fstat[bin] <- model$fstatistic[1]
      abs_beta[bin] <- abs(model$coefficients[2, 1])
      lm_pval[bin] <- model$coefficients[2,4]
      
    }
    
    ## Find out the longest stretch in which pval < 0.05 and create sum-of-F-stats statistic ##
    stretch <- lm_pval < 0.05
    indices <- stretch_start_end(stretch)  # what does stretch_start_end do???
    
    # if there is no stretch, take the max beta and max f statistic #
    if(is.na(indices[1])) {
      beta_stretch <- max(abs_beta)
      fstat_stretch <- max(fstat)
      # if there is a stretch, sum the betas and the f statistics #
    } else {
      beta_stretch <- sum(abs_beta[indices[1]:indices[2]]) # Summary stat
      fstat_stretch <- sum(fstat[indices[1]:indices[2]]) # Summary stat
    }
    
    
    # Create null distribution by shuffling labels #
    null_fstat_stretch <- NULL
    null_fstat_stretch <- foreach(h = 1:niter, .inorder=FALSE) %dopar% {
      # print status #
      if(h %% 100 == 0){
        print(paste0( (h/niter) * 100, "% complete for electrode", which(electrodes %in% elec), " of ", length(electrodes)))
      }
      
      # intialize vars #
      null_fstat <- NULL
      null_abs_beta <- NULL
      null_lm_pval <- NULL
      
     for (bin in nBins) { # This is the slooooow step
        
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
      
      return(null_fstat_stretch[h])
    }
    
    # our test statistic for the permutation test is the sum of f tests #
    perm_pval <- sum(null_fstat_stretch > fstat_stretch)/niter

    # save vals in results df #
    results <- results %>% 
      mutate_cond(electrode == elec, R2 = r2, absBeta = abs_beta, Fstat = fstat, p = lm_pval, fstretch = fstat_stretch, perm_p = perm_pval)
    
  
    # save results to results folder #
    write.csv(results, path(here(), "results", sub, paste0(region_name, "_", elec, "_", regressor, "_", tag, "_results.csv")))
  
  }
  
}
