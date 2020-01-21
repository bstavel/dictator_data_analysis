compile_results <- function(regions){
# function rbinds result dfs based on list of regions #
  
  # search the results folder for results csvs #
  filenames <- dir(fs::path(here(), "results"))
  files_to_combine <- filenames[grepl(paste(paste0("^", regions, "_.*results.csv"), collapse = "|"), filenames, perl = F)]
  
  # loop over and rbind #
  combined_results <- NULL
  for(file in files_to_combine){
    temp <- read.csv(fs::path(here(), "results", file))
    predictor <- gsub("^[^_]*_", "", gsub(".*-", "", gsub("_results.csv", "", file))) # because _someone_ didn't save the predictor in the damn df
    temp$predictor <- predictor
    combined_results <- rbind(combined_results, temp)
  }
  
  return(combined_results)
}