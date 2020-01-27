stepwise_anova_method_par <- function(regression_results, hg_behave, niter = 10000){
  
  # only active electrodes #
  active_results <- regression_results %>% filter(perm_p < (0.05)/length(regression_results$electrodes))
  
  # gets vars from df #
  electrodes <- unique(active_results$electrode)
  predictors <- unique(active_results$predictor)
  bins <- unique(regression_results$bin)
  
  
  for(elec in electrodes){
    # print status # 
    print(paste0("Beginning stepped anova analysis for electrode ", which(electrodes %in% elec), " out of ", length(electrodes)))
    
    # subset to the regressions for given elec #
    elec_results <- active_results %>% filter(electrode == elec)
    hg_elec <- hg_behave %>%filter(as.character(electrodes) == as.character(elec))
    
    # select predictor with highest R2 is this really what Ignacio does??? #
    max_r2_df <- elec_results %>% 
      group_by(predictor) %>%
      mutate(max_r2 = max(R2))
    
    # find max of maxes #
    chosen_predictors <- unique(max_r2_df$predictor[which.max(max_r2_df$max_r2)])
    
    # pool other predictors #
    remaining_predics <- predictors[!predictors %in% chosen_predictors]
    while(length(remaining_predics) > 1) { # have to figure this out
      f_anova <- matrix(NA, nrow = length(remaining_predics), ncol = length(bins))
      rownames(f_anova) <- remaining_predics
      colnames(f_anova) <- bins
      for(bin in bins) {
        # create df of only chosen regressors
        base_df <- hg_behave %>% select(bin, chosen_predictors)
        # change name so we don't have to paste0 #
        colnames(base_df)[grepl("bin", colnames(base_df))] <- "bin" # a bit neurotic but okay
        # run base model #
        base_model <- lm(bin ~ ., data = base_df)

        # iterate over other complex models #
        for(pred in remaining_predics){
          # create df of only chosen regressors + pred
          stepped_df <- hg_behave %>% select(bin, chosen_predictors, pred)
          
          # change name so we don't have to paste0 #
          colnames(stepped_df)[grepl("bin", colnames(stepped_df))] <- "bin" # a bit neurotic but okay
          
          # run base model #
          stepped_model <- lm(bin ~ ., data = stepped_df)

          # run anova #
          anova_temp <- anova(base_model, stepped_model)

          # store F val #
          f_anova[pred, bin] <- anova_temp$F[2]
        }
      }

      # get sig level via df of model #
      fstat_thr = qf(0.95, df1=anova_temp$Df[2], df2=anova_temp$Res.Df[2])

      # get sum of f stat for all predictors #
      fstat_stretch_temp <- NULL
      for(pred in remaining_predics) {
        stretch <- f_anova[pred,] > fstat_thr
        indices <- stretch_start_end(stretch)
        if(is.na(indices[1])) {
          fstat_stretch_temp[pred] <- max(f_anova[pred,])
        } else {
          fstat_stretch_temp[pred] <- sum(f_anova[pred,indices[1]:indices[2]]) # Summary stat
        }
      }

      # pick next predictor based on sum of f stat
      selected_unchosen_reg <-remaining_predics[which(fstat_stretch_temp == max(fstat_stretch_temp,na.rm=T))]
      chosen_predictors <- c(chosen_predictors, selected_unchosen_reg)
      remaining_predics <- remaining_predics[!(remaining_predics %in% chosen_predictors)]
      } # end while loop

    
    # Now run model comparisons in order to get the statistics #
    base_pred <- chosen_predictors[1]
    stepped_pred <- chosen_predictors[!chosen_predictors %in% base_pred]
    fstat_stretch <- NULL
    fstat_stretch_shuffle <- matrix(NA, nrow = length(stepped_pred), ncol = niter)
    rownames(fstat_stretch_shuffle) <- stepped_pred
    colnames(fstat_stretch_shuffle) <- 1:niter
    perm_p_anova <- NULL
    for(pred in stepped_pred) {
      f_anova <- NULL
      for(bin in bins){
      ## base model ##
      # create df of only chosen regressors
      base_df <- hg_behave %>% select(bin, base_pred)
      # change name so we don't have to paste0 #
      colnames(base_df)[grepl("bin", colnames(base_df))] <- "bin" # a bit neurotic but okay
      # run base model #
      base_model <- lm(bin ~ ., data = base_df)
      
      ## stepped model ##
      # create df of only chosen regressors + pred
      stepped_df <- hg_behave %>% select(bin, base_pred, pred)
      # change name so we don't have to paste0 #
      colnames(stepped_df)[grepl("bin", colnames(stepped_df))] <- "bin" # a bit neurotic but okay
      # run base model #
      stepped_model <- lm(bin ~ ., data = stepped_df)  
      
      ## anova ##
      anova_temp <- anova(base_model, stepped_model)
      f_anova[bin] <- anova_temp$F[2]
      }
      
      # get summed f statistic #
      # get sig level via df of model #
      fstat_thr = qf(0.95, df1=anova_temp$Df[2], df2=anova_temp$Res.Df[2])
      # get sum of f stat  #
        stretch <- f_anova > fstat_thr
        indices <- stretch_start_end(stretch)
        if(is.na(indices[1])) {
          fstat_stretch[pred] <- max(f_anova)
        } else {
          fstat_stretch[pred] <- sum(f_anova[indices[1]:indices[2]]) # Summary stat
        }
        
      ## run null models ##
        fstat_stretch_shuffle[pred, ] <- foreach(h = 1:niter, .inorder=FALSE, .combine = 'c') %dopar% {
        f_anova_shuffle <- NULL
        if(h %% 1 == 0){
          print(paste0( (h/niter) * 100, "% complete for predictor ", pred, ", on electrode ", which(electrodes %in% elec), " of ", length(electrodes)))
        }
        for(bin in bins){
          ## base model ##
          set.seed(h) # reproducability for anovas
          # create df of only chosen regressors and shuffle predictors
          base_df <- hg_behave %>% select(bin, base_pred) %>% mutate_at(vars(base_pred), sample)
          # change name so we don't have to paste0 #
          colnames(base_df)[grepl("bin", colnames(base_df))] <- "bin" # a bit neurotic but okay
          # run base model #
          base_model <- lm(bin ~ ., data = base_df)
          
          ## stepped model ##
          set.seed(h) # reproducability for anovas
          # create df of only chosen regressors + pred
          stepped_df <- hg_behave %>% select(bin, base_pred, pred) %>% mutate_at(vars(base_pred, pred), sample)
          # change name so we don't have to paste0 #
          colnames(stepped_df)[grepl("bin", colnames(stepped_df))] <- "bin" # a bit neurotic but okay
          # run base model #
          stepped_model <- lm(bin ~ ., data = stepped_df)  
          
          ## anova ##
          anova_temp <- anova(base_model, stepped_model)
          f_anova_shuffle[bin] <- anova_temp$F[2]
        } # end bins loop
        
        # get summed f statistic #
        # get sig level via df of model #
        fstat_thr = qf(0.95, df1=anova_temp$Df[2], df2=anova_temp$Res.Df[2])
        # get sum of f stat  #
        stretch <- f_anova_shuffle > fstat_thr
        indices <- stretch_start_end(stretch)
        if(is.na(indices[1])) {
          fstat_stretch_shuffle[h] <- max(f_anova_shuffle)
        } else {
          fstat_stretch_shuffle[h] <- sum(f_anova_shuffle[indices[1]:indices[2]]) # Summary stat
        }
        
        return(fstat_stretch_shuffle[h])
      } # end permutation
      
      # get permuted p value #
      perm_p_anova[pred] <- sum(fstat_stretch[pred] < fstat_stretch_shuffle[pred, ])/niter
      
      # update base predictors #
      base_pred <- c(base_pred, pred)
      
    } # end predictor loop
      
    # save results #
    stepped_anova_results <- data.frame(matrix(NA, nrow = (length(chosen_predictors) - 1), ncol = 5))
    colnames(stepped_anova_results) <- c("region", "electrode", "first_pred", "stepped_pred", "anova_perm_p")
    region_name <- unique(regression_results$region[regression_results$electrode == elec])
    stepped_anova_results$region <- region_name
    stepped_anova_results$electrode <- elec
    stepped_anova_results$first_pred <- chosen_predictors[1]
    stepped_anova_results$stepped_pred <- names(perm_p_anova)
    stepped_anova_results$anova_perm_p <- perm_p_anova
    write.csv(stepped_anova_results, path(here(), "results", paste0(region_name, "_", elec, "_stepped_anova_results.csv")))
    
    } # elec loop
    
} # function end
  