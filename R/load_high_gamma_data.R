load_high_gamma_data <- function(file_path_to_hg_data, file_path_to_electrode_names) {
  # load data #
  hg_data <- read.csv(file_path_to_hg_data, stringsAsFactors = F, header = F, comment.char="") # comment.char for fast loading
  # rename columns #
  hg_data <- hg_data %>%
    rename(index = V1750) %>%
    rename(trial = V1751)
  
  colnames(hg_data) <- gsub("^V", "time_", colnames(hg_data))
  
  # merge with electrode data #
  electrode_names <- read.csv(file_path_to_electrode_names)
  hg_data <- merge.data.frame(hg_data, electrode_names, by = "index", all = T)
  hg_data <- hg_data %>% select(-index) %>% rename(electrodes = Var1)
  
  return(hg_data)
}