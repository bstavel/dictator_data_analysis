load_high_gamma_data <- function(file_path_to_hg_data, file_path_to_electrode_names) {
  # load data #
  hg_data <- read.csv(file_path_to_hg_data, stringsAsFactors = F, header = F, comment.char="") # comment.char for fast loading
  # rename columns #
  index2 <- colnames(hg_data)[ncol(hg_data)]
  index1 <- colnames(hg_data)[ncol(hg_data)-1]
  hg_data <- hg_data %>%
    rename(index = all_of(index1)) %>%
    rename(trial = all_of(index2))
  
  colnames(hg_data) <- gsub("^V", "time_", colnames(hg_data))
  
  # merge with electrode data #
  electrode_names <- read.csv(file_path_to_electrode_names)
  hg_data <- merge.data.frame(hg_data, electrode_names, by = "index", all = T)
  hg_data <- hg_data %>% select(-index) %>% rename(electrodes = Var1)
  
  return(hg_data)
}

clean_electrode_names <- function(df){
  # clean the names
  df <- df %>%
    mutate(electrodes = gsub("POL ", "", electrodes)) %>%
    mutate(electrodes = gsub(" POL", "", electrodes)) %>%
    mutate(electrodes = gsub("-Ref", "", electrodes)) %>%
    mutate(electrodes = gsub("-Ref-", "-", electrodes)) 
  
  return(df)
}

load_high_gamma_data_clean <- function(file_path_to_hg_data, 
                                 file_path_to_electrode_names,
                                 file_path_to_elecs_of_interest) {
  # load data #
  hg_data <- read.csv(file_path_to_hg_data, stringsAsFactors = F, header = F, comment.char="") # comment.char for fast loading
  # rename columns #
  index2 <- colnames(hg_data)[ncol(hg_data)]
  index1 <- colnames(hg_data)[ncol(hg_data)-1]
  hg_data <- hg_data %>%
    rename(index = all_of(index1)) %>%
    rename(trial = all_of(index2))
  
  colnames(hg_data) <- gsub("^V", "time_", colnames(hg_data))
  
  # merge with electrode data #
  electrode_names <- read.csv(file_path_to_electrode_names)
  hg_data <- merge.data.frame(hg_data, electrode_names, by = "index", all = T)
  hg_data <- hg_data %>% select(-index) %>% rename(electrodes = Var1)
  
  # read in elecs #
  elecs_to_use <- read.csv(file_path_to_elecs_of_interest)
  
  # filter to good electrodes #
  hg_data <- hg_data %>% filter(grepl(paste(elecs_to_use$Electrode, collapse = "|"), electrodes))
  
  # clean names #
  hg_data <- clean_electrode_names(hg_data)
  
  return(hg_data)
}
