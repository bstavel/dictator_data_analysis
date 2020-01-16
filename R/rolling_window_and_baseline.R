rolling_window_and_baseline <- function(df, lWin = 100, lOver = 50){
  ### function to calculate the the mean in window specified as lWin with overlap lOver. ###
  #### then subtract the baseline, which is the first window. ###
  
  # First separate out indices #
  indices <- df %>% select(-starts_with("bin"))
  hg_df <- df %>% select(starts_with("bin"))

  # calculate the rolling average
  df_rollmean <- apply(hg_df, 1, function(x) rollapply(x, lWin, mean, by = lOver, align = "right", partial = F, by.column = T))
  df_rollmean <- data.frame(t(df_rollmean))
  
  # rename columns #
  colnames(df_rollmean) <- gsub("^X", "bin_", colnames(df_rollmean))
  
  # get baseline #
  baseline <- df_rollmean %>% mutate(baseline = (bin_1 + bin_2)/2) %>% select(baseline)
  
  # subtract the baseline (time around beginning of option presentation)
  df_rollmean_baseline <- apply(df_rollmean, 2, function(col) col - baseline$baseline)
  
  # rebind back #
  hg_clean <- cbind(indices, df_rollmean_baseline)
  
  # delete cols before time zero #
  hg_clean <- hg_clean %>% select(-bin_1, -bin_2)
  
 return(hg_clean)
}