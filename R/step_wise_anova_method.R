stepwise_anova_method <- function(regression_results, hg_behave){
  
  # only active electrodes #
  active_results <- regression_results %>% filter(perm_p < (0.05)/length(regression_results$electrodes))
  
  # gets vars from df #
  electrodes <- unique(active_results$electrode)
  predictors <- unique(active_results$predictor)
  bins <- unique(active_results$bin)
  
  
  for(elec in electrodes){
    
    if(length(predictors) < 2) {
      
      print(paste0("Zero or ! significant regressor for this electrode: ", elec))
      
    } else {
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
      }
      
    }
    
  }
  
}