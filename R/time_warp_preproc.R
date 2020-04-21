time_warp_preproc <- function(pres_data, choice_data, behave_data){
  
  # merge presentation data with behavioral data, and the select RT, trial, electrodes, and HFB data #
  pres_behave <- merge.data.frame(pres_data, behave_data, by.x = "trial", by.y = "round", all.x = T, all.y = F)
  pres_behave <- pres_behave %>% select(starts_with("time"), RT, trial, electrodes)
  
  # cut into three groups, baseline, pre-choice data, and post choice data #
  hg_behave_trial_parts <- pres_behave %>%
    gather(key = "time", value = "gamma", -RT, -trial, -electrodes) %>%
    mutate(time = as.numeric(gsub("time_", "", time))) %>%
    group_by(trial) %>%
    mutate(RT_cut = if_else(time < 200, "baseline", if_else(time < ( (RT*1000) + 200), "during", "post")))
  
  # average baseline data #
  baseline_averaged <- hg_behave_trial_parts %>%
    filter(RT_cut == "baseline") %>%
    group_by(trial, electrodes) %>%
    summarise_at(vars(gamma), mean) %>%
    rename(bin_0 = gamma) 
  
  # take the pre-choice data, and cut into 15 groups. Average the HFB data within groups #
  averaged_trial_parts <- hg_behave_trial_parts %>%
    filter(RT_cut == "during") %>%
    mutate(length_of_bins = ceiling((RT*1000)/15)) %>% # must be ceiling for RTs < 1s
    group_by(trial, electrodes) %>%
    mutate(averaging_groups_size  = n()) %>%
    mutate(averaging_groups = as.numeric(sapply(1:15, function(x) 
      rep(x, unique(length_of_bins))))[1:unique(averaging_groups_size)]) %>% # the 15th group will be smaller if RT doesn't divide into 15 evenly
    group_by(trial, electrodes, averaging_groups) %>%
    summarise_at(vars(gamma), mean) 
  
  # divide by baseline, unsure about this step given the bins size is variable #
  to_baseline_trial_bins <- merge.data.frame(averaged_trial_parts, baseline_averaged, by = c("trial", "electrodes"))
  baselined_trial_bins <- to_baseline_trial_bins %>%
    mutate(gamma = gamma - bin_0)
  
  # expand out #
  cleaned_pre_choice_data <- baselined_trial_bins %>%
    rename(bins = averaging_groups) %>%
    mutate(bins = paste0("bin_", bins)) %>%
    mutate(gamma = na_if(gamma, gamma < 2 & gamma > -2)) %>% 
    pivot_wider(id_cols = c(electrodes, trial, bin_0), 
                names_from = bins, 
                values_from = gamma)
  
  ## choice data ##
   
  # Separate out indices #
  indices <- choice_data %>% select(-starts_with("time"))
  hg_choice <- choice_data %>% select(starts_with("time"))
  hg_choice <- choice_data[, 750:2249] # exclude pre chocie time
  

  # calculate the rolling average, should this be left -- yes it is about the indexing anf the left side already has na vals so it is preferable to use that side
  lWin <- 100
  lOver <- 50
  df_rollmean <- apply(hg_choice, 1, function(x) rollapply(x, lWin, function(x) mean(x, na.rm = T), by = lOver, align = "left", partial = F, by.column = T))
  df_rollmean <- data.frame(t(df_rollmean))
  
  # rename, rebind back and expand #
  colnames(df_rollmean) <- paste0("bin_", 16:(ncol(df_rollmean) + 15))
  choice_meaned <- cbind(indices, df_rollmean)
  choice_meaned <-  choice_meaned %>%
    gather(key = "bins", value = "gamma", -trial, -electrodes)
  
  # subtract the baseline (time around beginning of option presentation)
  to_baseline_choice <- merge.data.frame(choice_meaned, baseline_averaged, by = c("trial", "electrodes"))
  baselined_choice_bins <- to_baseline_choice %>%
    mutate(gamma = gamma - bin_0)
  
  # expand out #
  cleaned_post_choice_data <- baselined_choice_bins %>%
    mutate(gamma = na_if(gamma, gamma < 2 & gamma > -2)) %>%
    pivot_wider(id_cols = c(electrodes, trial, bin_0), 
                names_from = bins, 
                values_from = gamma) %>%
    select(-bin_0)
    
  ## clean and add behave data ##
  # combine presentation and choice #
  clean_hfb_data <- merge.data.frame(cleaned_pre_choice_data, cleaned_post_choice_data, by = c("trial", "electrodes"))
  # merge with behavioral data #
  clean_hfb_behave_data <- merge.data.frame(clean_hfb_data, behave_data, by.x = "trial", by.y = "round", all.x = T, all.y = F)
  # clean any missing RTs #
  clean_hfb_behave_data <- clean_hfb_behave_data %>% filter(!is.na(RT))
  
  return(clean_hfb_behave_data)
  
}


