create_choice_time_data <- function(df, behave_data){
  
  # separate from indices #
  indicies <- df %>% select(trial, electrodes)
  hg_data <- df %>% select(starts_with("time"))
  # pull baseline #
  baseline_df <- hg_data[, 1:200]
  baseline <- apply(baseline_df, 1, function(x) mean(x, na.rm = T))
  # recombine #
  hg_data_tmp <- hg_data
  hg_data_baselined <- apply(hg_data_tmp, 2, function(col) t(as.vector(col - baseline)))
  hg_data_baselined <- as.data.frame(hg_data_baselined)
  hg_data_baselined <- cbind(hg_data_baselined, indicies)
  # merge with behavior #
  hg_raw_plot <- merge.data.frame(hg_data_baselined, behave_data, by.x = "trial", by.y = "round", all.x = T, all.y = F)
  
  # convert to choice locked #
  df_choice_locked <- hg_raw_plot %>%
    pivot_longer(names_to = "time", values_to = "gamma", cols = starts_with("time_")) %>%
    mutate(time = (as.numeric(gsub("time_", "", time))) -200 ) %>%
    group_by(electrodes, RT, trial) %>%
    mutate(time = time - RT*1000) %>%
    mutate(time = paste0("time_", time)) %>%
    pivot_wider(names_from = time, values_from = gamma)

  
}