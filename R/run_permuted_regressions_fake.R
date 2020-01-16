run_permuted_regressions <- function(brain_behave_data, electrodes, nBins, region_name, niter = 1000) {
  
  # create results dfs #
  results_advantageous <- data.frame(matrix(nrow = length(electrodes)*length(nBins), ncol = 8))
  results_disadvantageous <- data.frame(matrix(nrow = length(electrodes)*length(nBins), ncol = 8))
  
  colnames(results_advantageous) <- c("electrode", "bin", "R2", "absBeta", "Fstat", "p", "fstretch", "perm_p")
  colnames(results_disadvantageous) <- c("electrode", "bin", "R2", "absBeta", "Fstat", "p", "fstretch", "perm_p")
  
  # fill results dfs with electrode & bin info #
  results_advantageous$electrode <- sort(rep(as.character(electrodes), length(nBins)))
  results_advantageous$bin <- rep(nBins, length(electrodes))
  
  results_disadvantageous$electrode <- sort(rep(as.character(electrodes), length(nBins)))
  results_disadvantageous$bin <- rep(nBins, length(electrodes))
  
  
  for (elec in electrodes) {
    
    # filter to single electrode #
    brain_behave_data_elec <- brain_behave_data %>% filter(electrodes == elec)
    # brain_behave_data_elec <- brain_behave_data_elec[complete.cases(brain_behave_data_elec), ]
    
    # initialize temp vars #
    r2_disadv <- NULL
    fstat_disadv <- NULL
    abs_beta_disadv <- NULL
    lm_pval_disadv <- NULL
    r2_adv <- NULL
    fstat_adv <- NULL
    abs_beta_adv <- NULL
    lm_pval_adv <- NULL
    
    for (bin in nBins) {
      # run models #
      bin_vec <- brain_behave_data_elec[, bin]
      model_disadvant <- summary(lm(bin_vec ~ sample(1:228), data = brain_behave_data_elec))
      model_advant <- summary(lm(bin_vec ~ sample(1:228), data = brain_behave_data_elec))
      
      # store info from models #
      r2_disadv[bin] <- model_disadvant$r.squared
      fstat_disadv[bin] <- model_disadvant$fstatistic[1]
      abs_beta_disadv[bin] <- abs(model_disadvant$coefficients[2, 1])
      lm_pval_disadv[bin] <- model_disadvant$coefficients[2,4]
      
      r2_adv[bin] <- model_advant$r.squared
      fstat_adv[bin] <- model_advant$fstatistic[1]
      abs_beta_adv[bin] <- abs(model_advant$coefficients[2,1])
      lm_pval_adv[bin] <- model_advant$coefficients[2,4]
    }
    
    ## Find out the longest stretch in which pval < 0.05 and create sum-of-F-stats statistic ##
    # disadventageous #
    stretch <- lm_pval_disadv < 0.05
    indices <- stretch_start_end(stretch)  # what does stretch_start_end do???
    
    # if there is no stretch, take the max beta and max f statistic #
    if(is.na(indices[1])) {
      beta_stretch_disadv <- max(abs_beta_disadv)
      fstat_stretch_disadv <- max(fstat_disadv)
      # if there is a stretch, sum the betas and the f statistics #
    } else {
      beta_stretch_disadv <- sum(abs_beta_disadv[indices[1]:indices[2]]) # Summary stat
      fstat_stretch_disadv <- sum(fstat_disadv[indices[1]:indices[2]]) # Summary stat
    }
    
    # adventageous #
    stretch <- lm_pval_adv < 0.05
    indices <- stretch_start_end(stretch)  # what does stretch_start_end do???
    
    # if there is no stretch, take the max beta and max f statistic #
    if(is.na(indices[1])) {
      beta_stretch_adv <- max(abs_beta_adv)
      fstat_stretch_adv <- max(fstat_adv)
      # if there is a stretch, sum the betas and the f statistics #
    } else {
      beta_stretch_adv <- sum(abs_beta_adv[indices[1]:indices[2]]) # Summary stat
      fstat_stretch_adv <- sum(fstat_adv[indices[1]:indices[2]]) # Summary stat    
    }
    
    # Create null distribution by shuffling labels #
    null_beta_stretch_disadv <- NULL
    null_fstat_stretch_disadv <- NULL
    null_beta_stretch_adv <- NULL
    null_fstat_stretch_adv <- NULL
    for (h in 1:niter) {
      # intialize vars #
      null_fstat_disadv <- NULL
      null_abs_beta_disadv <- NULL
      null_lm_pval_disadv <- NULL
      
      null_fstat_adv <- NULL
      null_abs_beta_adv <- NULL
      null_lm_pval_adv <- NULL
      for (bin in nBins) { # This is the slooooow step
        
        # if(remove_outliers) {
        #   temp <- data.frame(peak=peak)
        #   temp <- remove_outliers(temp,peak)
        #   peak <- temp$peak
        # }
        
        # run models #
        bin_vec <- brain_behave_data_elec[, bin]
        null_lm_disadv <-summary(lm(bin_vec ~ sample(1:228), data = brain_behave_data_elec))
        null_lm_adv <- summary(lm(bin_vec ~ sample(1:228), data = brain_behave_data_elec))
        
        # save info from models #
        null_fstat_disadv[bin] <- null_lm_disadv$fstatistic[1]
        null_abs_beta_disadv[bin] <- abs(null_lm_disadv$coefficients[2,1])
        null_lm_pval_disadv[bin] <- null_lm_disadv$coefficients[2,4]
        
        null_fstat_adv[bin] <- null_lm_adv$fstatistic[1]
        null_abs_beta_adv[bin] <- abs(null_lm_adv$coefficients[2,1])
        null_lm_pval_adv[bin] <- null_lm_adv$coefficients[2,4]      
      }
      
      ## Find out the longest stretch in which pval < 0.05 and create sum-of-F-stats statistic ##
      # disadvantageous #
      stretch <- null_lm_pval_disadv < 0.05
      indices <- stretch_start_end(stretch)
      if(is.na(indices[1])) {
        # if no stretch just take the max vals #
        null_beta_stretch_disadv[h] <- max(null_abs_beta_disadv)
        null_fstat_stretch_disadv[h] <- max(null_fstat_disadv)
        
      } else {
        # if a stretch take the sum of f stats #  
        null_beta_stretch_disadv[h] <- sum(null_abs_beta_disadv[indices[1]:indices[2]]) # Summary stat
        null_fstat_stretch_disadv[h] <- sum(null_fstat_disadv[indices[1]:indices[2]]) # Summary stat
      }
      
      # advantageous #
      stretch <- null_lm_pval_adv < 0.05
      indices <- stretch_start_end(stretch)
      if(is.na(indices[1])) {
        # if no stretch just take the max vals #
        null_beta_stretch_adv[h] <- max(null_abs_beta_adv)
        null_fstat_stretch_adv[h] <- max(null_fstat_adv)
        
      } else {
        # if a stretch take the sum of f stats #  
        null_beta_stretch_adv[h] <- sum(null_abs_beta_adv[indices[1]:indices[2]]) # Summary stat
        null_fstat_stretch_adv[h] <- sum(null_fstat_adv[indices[1]:indices[2]]) # Summary stat
        
      }    
    }
    
    # our test statistic for the permutation test is the sum of f tests #
    perm_pval_disadv <- sum(null_fstat_stretch_disadv > fstat_stretch_disadv)/niter
    perm_pval_adv <- sum(null_fstat_stretch_adv > fstat_stretch_adv)/niter
    
    if(perm_pval_adv == 0 | perm_pval_disadv == 0) {
      stop("p value nonsense!")
    }
    
    # save vals in results df #
    results_advantageous <- results_advantageous %>% 
      mutate_cond(electrode == elec, R2 = r2_adv, absBeta = abs_beta_adv, Fstat = fstat_adv, p = lm_pval_adv, fstretch = fstat_stretch_adv, perm_p = perm_pval_adv)
    
    results_disadvantageous <- results_disadvantageous %>% 
      mutate_cond(electrode == elec, R2 = r2_disadv, absBeta = abs_beta_disadv, Fstat = fstat_disadv, p = lm_pval_disadv, fstretch = fstat_stretch_disadv, perm_p = perm_pval_disadv)
    
    
  }
  
  region_name <- "fake"
  # save results to results folder #
  write.csv(results_disadvantageous, path(here(), "results", paste0(region_name, "_results_disadvantageous.csv")))
  write.csv(results_advantageous, path(here(), "results", paste0(region_name, "_results_advantageous.csv")))
  
}
