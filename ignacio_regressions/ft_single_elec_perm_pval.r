ft_single_elec_perm_pval = function(SUBID, elec_name, reg, freq_band = 'hg', epoch = 'buttonpress', remove_outliers = F, niter = 1000) {
  
  # Generates a permutation Null distribution for a given electrode, regressor and frequency band
  # 1. Carry out the linear regression
  # 2. Obtain f-stat threshold for a=0.05
  # 3. Obtain max # consecutive datapoints with p < 0.05
  # 4. Summary statistic: sum of abs(betas) across all those datapoints
  # 5. Calculate same statistic over niter label permutations: Null
  # 6. Generate and return perm pval
  # Check cluster or local

  # Initialize
  wd = getwd()
  wd = substr(wd,1,10)
  path = '/Users/ignaciosaez/Dropbox/Berkeley/Project_EcoG'
  lpfc_path = '/Users/ignaciosaez/Dropbox/Berkeley/Project_EcoG_LPFC'
  library(R.matlab)
  source(paste(path,'/scripts/functions/lmp.r',sep=''))
  source(paste(path,'/scripts/functions/stretch_start_end.r',sep=''))
  source(paste(path,'/scripts/functions/ofc_to_subj.r',sep=''))
  source(paste(path,'/scripts/functions/get_reg_matrix.r',sep=''))
  source(paste(path,'/scripts/functions/remove_outliers_IQR.r',sep=''))
  freq_band = tolower(freq_band)
  rotate <- function(x) t(apply(x, 2, rev))
  SUBID_list = c('s01_GP47','s02_ST35','s03_AMC046','s04_GP48','s05_ST40','s06_IR06','s07_GP51','s08_IR09','s09_IR10','s10_IR14','s11_IR16','s12_IR17','s13_IR19','s14_IR20','s15_IR21','s16_IR22','s17_IR26','s18_ST46','s19_IR27','s20_IR28','s21_IR31','s22_IR35','s23_OS7','s24_IR36','s25_IR38','s26_OS11','s27_IR39','s28_IR40','s29_OS16','s30_IR41','s31_IR53','s32_IR55','s33_IR56')
  SUBID_behav = paste('gamble.data.s',sprintf('%02d', which(SUBID_list == SUBID)),'.csv',sep='')
  
  # Catch error
  if(is.null(reg)) {
    print('reg is empty - cannot run single_elec_perm_pval')
    return()
  }
  reg = gsub('_','.',reg)
  
  # Info message
  print(paste('Running single_elec_perm_pval_ft on ',SUBID,', elec ',elec_name,', reg ',reg,', ',epoch,' epoch, ',freq_band,' freq_band',sep=''))
  #print('Running lm sum-of-Fstat permutation statistics')
  if(remove_outliers) {
    print('Removing outliers')
  }
  
  # Load HG/behavioral data
  behav_globals = readMat(paste(path,'/data/',SUBID,'/behav_globals.mat',sep=''))
  
  # Read parsed windowed data
  events = readMat(paste(lpfc_path, '/',SUBID,'/parsed/ft_parsed_',epoch,'_',elec_name,'.mat',sep=''))
  events = events$window.events[,,1]
  eval(parse(text=paste('events = events$',tolower(freq_band),sep='')))
  
  # Get parameters
  nevents = dim(events)[1]
  nwin = dim(events)[2]
  reg_temp = eval(parse(text=paste('behav_globals$',reg,sep='')))
  z_reg_temp = scale(reg_temp)[,1]
  if(SUBID == 's32_IR55') { z_reg_temp = z_reg_temp[1:135] }
  
  # Run linear regression: get f-stat for each time window, fstat_thr and longest (fstat > fstat_thr) stretch
  rsq = fstat = abs_beta = lm_pval = NaN
  for (i in 1:nwin) {
    peak = events[,i]
    if(remove_outliers) {
      temp = data.frame(peak=peak)
      temp = remove_outliers_IQR(temp,peak)
      peak = temp$peak
    }
    temp_lm = summary(lm(peak ~ z_reg_temp))
    fstat[i] = temp_lm$fstatistic[1]
    abs_beta[i] = abs(temp_lm$coefficients[2,1])
    lm_pval[i] = temp_lm$coefficients[2,4]
  }
  
  # Find out the longest stretch in which pval > 0.05 and create sum-of-F-stats statistic
  stretch = lm_pval < 0.05
  indices = stretch_start_end(stretch)
  if(is.na(indices[1])) {
    beta_stretch = max(abs_beta)
    fstat_stretch = max(fstat)
  } else {
    beta_stretch = sum(abs_beta[indices[1]:indices[2]]) # Summary stat
    fstat_stretch = sum(fstat[indices[1]:indices[2]]) # Summary stat
  }
  
  # Same, after shuffling labels
  beta_stretch_null = fstat_stretch_null = fstat_null = 0
  for(h in 1:niter) {
    abs_beta_null = lm_pval_null = 0
    for (i in 1:nwin) { # This is the slooooow step
      peak = events[,i]
      if(remove_outliers) {
        temp = data.frame(peak=peak)
        temp = remove_outliers(temp,peak)
        peak = temp$peak
      }
      temp_lm = summary(lm(peak ~ sample(z_reg_temp)))
      fstat_null[i] = temp_lm$fstatistic[1]
      abs_beta_null[i] = abs(temp_lm$coefficients[2,1])
      lm_pval_null[i] = temp_lm$coefficients[2,4]
    }
    
    # Find out the longest stretch in which pval > 0.05 and create sum-of-F-stats statistic
    stretch = lm_pval_null < 0.05
    indices = stretch_start_end(stretch)
    if(is.na(indices[1])) {
      beta_stretch_null[h] = max(abs_beta_null)
      fstat_stretch_null[h] = max(fstat_null)
    } else {
      beta_stretch_null[h] = sum(abs_beta_null[indices[1]:indices[2]]) # Summary stat
      fstat_stretch_null[h] = sum(fstat_null[indices[1]:indices[2]]) # Summary stat
    }
    # if(fstat_stretch_null[h] > 140) { break }
  }
  
  #perm_pval = sum(beta_stretch_null > beta_stretch)/niter
  perm_pval = sum(fstat_stretch_null > fstat_stretch)/niter
  
  return(perm_pval)
  
}