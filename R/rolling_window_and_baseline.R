rolling_window_and_baseline <- function(df, lWin = 100, lOver = 50, choice_locked){
  ### function to calculate the the mean in window specified as lWin with overlap lOver. ###
  #### then subtract the baseline, which is the first window. ###
  
  if(choice_locked == TRUE) {
  
  # First separate out indices #
  indices <- df %>% select(-starts_with("time"))
  hg_df <- df %>% select(starts_with("time"))
  
  # then separate out pre trial baseline, taken from the presentation locked trials #
  baseline_df <- read.csv(path(here(), "munge", "baseline_data.csv"))
  baseline_df <- baseline_df %>% select(-X)
  
  # # reorder nas #
  # right_shift_na <- function(x) {
  #   na_indices <- which(is.na(x))
  #   not_na <-  which(!is.na(x))
  #   if(all(min(na_indices) > not_na)){ # all of the na indicies should be greater than the non na indices
  #     new_x <- c(x[na_indices], x[not_na])
  #   } else {
  #     print(x[1])
  #   }
  #   return(unlist(new_x))  # will error if it is not pass the if all test
  # }
  # df_reorder <- apply(hg_df, 1, right_shift_na)
  # df_reorder <- data.frame(t(df_reorder))

  # calculate the rolling average, should this be left -- yes it is about the indexing anf the left side already has na vals so it is preferable to use that side
  df_rollmean <- apply(hg_df, 1, function(x) rollapply(x, lWin, function(x) mean(x, na.rm = T), by = lOver, align = "left", partial = F, by.column = T))
  df_rollmean <- data.frame(t(df_rollmean))

  # calculate rolling average for baseline 
  baseline <- apply(baseline_df, 1, function(x) mean(x, na.rm = T))
  baseline <- data.frame(t(baseline))
  
  # rename columns #
  colnames(df_rollmean) <- c(rev(paste0("pre_", 1:13)),  paste0("post_", 1:20))
  
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
    write.csv(baseline, path(here(), "munge", "baseline.csv"))

    hg_df <- hg_df %>% select(201:ncol(hg_df))
    
    # calculate the rolling average
    df_rollmean <- apply(hg_df, 1, function(x) rollapply(x, lWin, mean, by = lOver, align = "left", partial = F, by.column = T))
    df_rollmean <- data.frame(t(df_rollmean))
    
    # rename columns #
    colnames(df_rollmean) <- gsub("^X", "bin_", colnames(df_rollmean))
    
    # subtract the baseline (time around beginning of option presentation)
    df_rollmean_baseline <- apply(df_rollmean, 2, function(col) t(as.vector(col - baseline)))
    
    # rebind back #
    hg_clean <- cbind(indices, df_rollmean_baseline)
    
    # delete cols before time zero #
    hg_clean <- hg_clean %>% select(-bin_1, -bin_2)
    
  }
 
  return(hg_clean) 
  
}  
