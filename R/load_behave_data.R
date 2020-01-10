load_behave_data = function(datapath){
  filenames <- dir(datapath,"behav.csv")
  DATA <- NULL
  
  for (i in 1:length(filenames)) {
    temp <- cbind(SID=substr(filenames[i],1,6),read.csv(fs::path(datapath, filenames[i]), header=TRUE))
    temp$RT <- as.numeric(as.character(temp$RT))
    DATA <- rbind(DATA,temp)
  }
  
  write.csv(DATA, fs::path(here(), "munge", "combined_behavioral_data.csv"))
  
}