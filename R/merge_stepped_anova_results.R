merge_stepped_anova_results <- function(regions){
 # function rbinds result dfs based on list of regions #
  
  # search the results folder for results csvs #
  filenames <- dir(fs::path(here(), "results"))
  files_to_combine <- filenames[grepl(paste(paste0("^", regions, "_.*results.csv"), collapse = "|"), filenames, perl = F)]
  files_to_combine <- files_to_combine[grepl("stepped_anova", files_to_combine)] # HACK fix this later please
  
  # loop over and rbind #
  combined_results <- NULL
  for(file in files_to_combine){
    temp <- read.csv(fs::path(here(), "results", file))
    combined_results <- rbind(combined_results, temp)
  }
  
  return(combined_results)
  
  
}