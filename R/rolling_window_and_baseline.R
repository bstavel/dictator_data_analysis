rolling_window_and_baseline <- function(df, baseline_csv, lWin = 150, lOver = 50, choice_locked){
  ### function to calculate the the mean in window specified as lWin with overlap lOver. ###
  #### then subtract the baseline, which is the first window. ###
  
  if(choice_locked == TRUE) {
    
  # order electrodes and trials to ensure good match with baseline #
  df <- df %>% arrange(electrodes, trial)
  
  # Separate out indices #
  indices <- df %>% select(-starts_with("time"))
  hg_df <- df %>% select(starts_with("time"))
  
  # then separate out pre trial baseline, taken from the presentation locked trials #
  baseline_df <- read.csv(baseline_csv)
  baseline_df <- baseline_df %>% arrange(hg_data.electrodes, hg_data.trial)
  baseline_mean <- baseline_df %>% select(starts_with("time")) %>% apply(., 1, function(x) mean(x))
  baseline <- data.frame(t(baseline_mean))

  # calculate the rolling average, should this be left -- yes it is about the indexing anf the left side already has na vals so it is preferable to use that side
  df_rollmean <- apply(hg_df, 1, function(x) rollapply(x, lWin, function(x) mean(x, na.rm = T), by = lOver, align = "left", partial = F, by.column = T))
  df_rollmean <- data.frame(t(df_rollmean))

  # rename columns #
  colnames(df_rollmean) <- c(rev(paste0("pre_", 1:13)),  paste0("post_", 1:(ncol(df_rollmean)-13)))
  # colnames(df_rollmean) <- c(rev(paste0("pre_", 1:5)),  paste0("post_", 1:10))
  
  # subtract the baseline (time around beginning of option presentation)
  df_rollmean_baseline <- apply(df_rollmean, 2, function(col) t(as.vector(col - baseline)))
  
  # rebind back #
  hg_clean <- cbind(indices, df_rollmean_baseline)
  
  } else if(choice_locked == FALSE) {
    
    ### function to calculate the the mean in window specified as lWin with overlap lOver. ###
    #### then subtract the baseline, which is the first window. ###
    
    # First separate out indices #
    indices <- df %>% select(-starts_with("time"))
    hg_df <- df %>% select(starts_with("time"))
    
    baseline_df <- hg_df %>% select(1:200)
    baseline <- apply(baseline_df, 1, function(x) mean(x, na.rm = T))
   
    hg_df <- hg_df %>% select(201:ncol(hg_df))
    
    # calculate the rolling average
    df_rollmean <- apply(hg_df, 1, function(x) rollapply(x, lWin, mean, by = lOver, align = "left", partial = F, by.column = T))
    df_rollmean <- data.frame(t(df_rollmean))
    
    # rename columns #
    colnames(df_rollmean) <- gsub("^X", "bin_", colnames(df_rollmean))
    
    # subtract the baseline (time around beginning of option presentation)
    df_rollmean_baseline <- apply(df_rollmean, 2, function(col) t(as.vector(col - baseline)))
    
    # rebind back # Note, we do not add back the baseline because it would be zero since we baseline with a 200 ms window AND average with a 200 second window
    hg_clean <- cbind(indices, df_rollmean_baseline)
    
  }
 
  return(hg_clean) 
  
}  
