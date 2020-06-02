load_all_freq_data <- function(file_path_to_tf_data, file_path_to_electrode_names) {
  # load data #
  tf_data <- read.csv(file_path_to_tf_data, stringsAsFactors = F, header = F, comment.char="") # comment.char for fast loading
  # rename columns #
  index2 <- colnames(tf_data)[ncol(tf_data) -1]
  index1 <- colnames(tf_data)[ncol(tf_data)-2]
  index3 <- colnames(tf_data)[ncol(tf_data)]
  tf_data <- tf_data %>%
    rename(index = index1) %>%
    rename(trial = index2) %>% 
    rename(frequency = index3)
  
  colnames(tf_data) <- gsub("^V", "time_", colnames(tf_data))
  
  # merge with electrode data #
  electrode_names <- read.csv(file_path_to_electrode_names)
  tf_data <- merge.data.frame(tf_data, electrode_names, by = "index", all = T)
  tf_data <- tf_data %>% select(-index) %>% rename(electrodes = Var1)
  
  return(tf_data)
}