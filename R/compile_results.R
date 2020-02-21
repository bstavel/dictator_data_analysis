compile_results <- function(regions, type){
# function rbinds result dfs based on list of regions #
  
  ## search the results folder for results csvs ##
  filenames <- dir(fs::path(here(), "results"))
  # only from the region list#
  files_to_combine <- filenames[grepl(paste(paste0("^", regions, "_.*results.csv"), collapse = "|"), filenames, perl = F)]
  # not anova results #
  files_to_combine <- files_to_combine[!grepl("anova", files_to_combine)] # HACK fix this later please
  # type #
  files_to_combine <- files_to_combine[grepl(type, files_to_combine)] # HACK fix this later please
  
  # loop over and rbind #
  combined_results <- NULL
  for(file in files_to_combine){
    temp <- read.csv(fs::path(here(), "results", file))
    predictor <- gsub("^[^_]*_", "", gsub(".*-", "", gsub(paste0("_", type, "_results.csv"), "", file))) # because _someone_ didn't save the predictor in the damn df
    region <- gsub("_.*", "", file)
    temp$predictor <- predictor
    temp$region <- region
    combined_results <- rbind(combined_results, temp)
  }
  
  return(combined_results)
}