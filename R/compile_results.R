compile_results <- function(regions, type, sub, path, spec_vars = c("*", "*"), absBeta = F){
# function rbinds result dfs based on list of regions #
  
  ## search the results folder for results csvs ##
  filenames <- dir(path)
  # only from the region list#
  files_to_combine <- filenames[grepl(paste(paste0("^", regions, "_.*results.csv"), collapse = "|"), filenames, perl = F)]
  # on certain variables? #
  files_to_combine <- files_to_combine[grepl(paste(spec_vars, collapse = "|"), files_to_combine)]
  # not anova results #
  files_to_combine <- files_to_combine[!grepl("anova", files_to_combine)] # HACK fix this later please
  # not beta results #
  files_to_combine <- files_to_combine[!grepl("full_beta", files_to_combine)] # HACK fix this later please
  # type #
  files_to_combine <- files_to_combine[grepl(type, files_to_combine)] # HACK fix this later please
  
  # loop over and rbind #
  combined_results <- NULL
  for(file in files_to_combine){
    temp <- read.csv(path(path, file))
    predictor <- gsub("^[^_]*_", "", gsub(".*-", "", gsub(paste0("_", type, "_results.csv"), "", file))) # because _someone_ didn't save the predictor in the damn df
    region <- gsub("_.*", "", file)
    temp$predictor <- predictor
    temp$region <- region
    # # grab not absolute value beta for exist check #
    not_abs_beta_file <- gsub("_results", "_full_beta_results", file)
    if(absBeta == T & file.exists(path(path, not_abs_beta_file))){
      # read in beta file #
      full_beta_tmp <- read.csv(path(path, not_abs_beta_file))
      # add beta col to temp df #
      temp$Beta <- full_beta_tmp$Beta
      temp <- temp %>% select(-absBeta)
    }
    # remove absolute beta so it can be combined with other values
    temp <- temp %>% select(-X)
    combined_results <- rbind(combined_results, temp)
  }
  
  return(combined_results)
}