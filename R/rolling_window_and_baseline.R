rolling_window_and_baseline <- function(df, lWin = 100, lOver = 50){
  ### function to calculate the the mean in window specified as lWin with overlap lOver. ###
  #### then subtract the baseline, which is the first window. ###
  
  # First separate out indices #
  indices <- df %>% select(-starts_with("time"))
  hg_df <- df %>% select(starts_with("time"))
  
  # then separate out pre trial baseline #
  baseline_df <- hg_df %>% select(1:199)
  task_df <- hg_df %>% select(200:ncol(hg_df))
  
  # reorder nas #
  right_shift_na <- function(x) {
    na_indices <- which(is.na(x))
    not_na <-  which(!is.na(x))
    if(all(min(na_indices) > not_na)){ # all of the na indicies should be greater than the non na indices
      new_x <- c(x[na_indices], x[not_na])
    } else {
      print(x[1])
    }
    return(unlist(new_x))  # will error if it is not pass the if all test
  }
  df_reorder <- apply(task_df, 1, right_shift_na)
  df_reorder <- data.frame(t(df_reorder))

  # calculate the rolling average, should this be left -- yes it is about the indexing anf the left side already has na vals so it is preferable to use that side
  df_rollmean <- apply(df_reorder, 1, function(x) rollapply(x, lWin, function(x) mean(x, na.rm = T), by = lOver, align = "left", partial = F, by.column = T))
  df_rollmean <- data.frame(t(df_rollmean))

  # calculate rolling average for baseline 
  baseline <- apply(baseline_df, 1, function(x) mean(x, na.rm = T))
  baseline <- data.frame(t(baseline))
  
  # rename columns #
  colnames(df_rollmean) <- c(rev(paste0("pre_", 1:40)),  paste0("post_", 1:20))
  
  # subtract the baseline (time around beginning of option presentation)
  df_rollmean_baseline <- apply(df_rollmean, 2, function(col) t(as.vector(col - baseline)))
  
  # rebind back #
  hg_clean <- cbind(indices, df_rollmean_baseline)
  
 return(hg_clean)
}