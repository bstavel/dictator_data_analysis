load_high_gamma_data <- function(file_path_to_hg_data, file_path_to_electrode_names) {
  # load data #
  hg_data <- read.csv(file_path_to_hg_data, stringsAsFactors = F, nrows = 18924, header = F, comment.char="")
  # rename columns #
  hg_data <- hg_data %>%
    rename(index = V3201) %>%
    rename(trial = V3202)
  colnames(hg_data) <- gsub("^V", "bin_", colnames(hg_data))
  
  # merge with electrode data #
  electrode_names <- read.csv(file_path_to_electrode_names)
  hg_data <- merge.data.frame(hg_data, electrode_names, by = "index", all = T)
  hg_data <- hg_data %>% select(-index) %>% rename(electrodes = Var1)
  
  return(hg_data)
}