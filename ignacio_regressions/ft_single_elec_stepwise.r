ft_single_elec_stepwise = function(SUBID, elec_name, reg_list, freq_band = 'hg', corr_nelecs = 200, epoch = 'outcome', pdf=T, niter=100, return_results=F) {
  # TO-DO list
  # 0 - FIX MULTIPLE COMPARISONS - now done as *#ofc_elecs in corr_pval, but the number of electrodes grows!!
  # 1. permutation pvals should be based on current round only - otherwise opening the time window = ++false positives!
  
  # After Padoa-Schioppa 2006: incrementally add regressors until no more variance is captured
  # Includes optimal regressor ordering, model selection and multiplexing.
  # Takes as input FT version of data (Nov 2018)
  
  # For each electrode:
  # 1. Run independent linear regressions for all regressors
  # 2. Select significant regressors:
  #   a. Restrict to any p<0.05
  #   b. Further restrict based on permutation test p<0.05
  #   c. Further restrict based on multiple comparison correction (by nelecs)
  # 3. Select the best regressor order:
  #   a. Start with regressor with highest RSq peak
  #   b. At each step, find and add component that adds most explanatory power (based on model comp ANOVA summary Fstat)
  # 4. Run stepwise regression: run permutation tests to select optimal number of regressors:
  #    Incrementally add regressors in optimal order until no more variance is captured (ANOVA p>0.05).
  # 5. Plot
  # 6. Save data and return as list
  
  # Permutation tests: used for regressor selection and model comparison.
  #   a. ANOVA Fstat timecourse is calculated (i.e. for simple-complex moddimel comparison)
  #   b. Summary Fstat is calculated = sum of Fstats over fstat_thr (corresponding to pval = 0.05)
  #   c. A Null is calculated by shuffling labels niter times
  #   d. Calculate and report permutation pvals

  # Fields in model output:
  # ncomp            Number of components in final model
  # comp_names       (1 x ncomp) Names of components in final model - ordered by best fit
  # comp_pos         (1 x length(reg_list)) Index of ordered components in original reg_list
  # lm_rsq_all       (nreg x time) Rsq for all regressors 
  # lm_pval_all      (nreg x time) pval for all regressors
  # lm_rsq           (nreg x time) RSq for individual linear regressions for kept regressors only
  # lm_pval          (nreg x time) pval for kept regressors
  # reg_pval_perm    corrected pvals for regressor selection (for all regressors that initially survived mult.comp.)
  # parse.params     parse parameters
  # if >2 components:
    # model_comp_pval_perm   (1 x ncomp) pvals for model comparison 
    # stepwise_model         (Rsq x time) Final model fit (incremental RSq)
    # residual_fits          (Rsq x time) Residual fits for stepwise regressors (residual RSq)
  
  ## 0. Initialize ####
  # Allow numeric SUBID
  if(is.numeric(SUBID)){
    SUBID_list = c('s01_GP47','s02_ST35','s03_AMC046','s04_GP48','s05_ST40','s06_IR06','s07_GP51','s08_IR09','s09_IR10','s10_IR14','s11_IR16','s12_IR17','s13_IR19','s14_IR20','s15_IR21','s16_IR22','s17_IR26','s18_ST46','s19_IR27','s20_IR28','s21_IR31','s22_IR35','s23_OS7','s24_IR36','s25_IR38','s26_OS11','s27_IR39','s28_IR40','s29_OS16','s30_IR41','s31_IR53','s32_IR55','s33_IR56')
    SUBID = SUBID_list[SUBID]
  }
  
  library(R.matlab)
  library(fields)
  library(RColorBrewer)
  
  # Check cluster or local
  path = '/Users/ignaciosaez/Dropbox/Berkeley/Project_EcoG'
  lpfc_path = '/Users/ignaciosaez/Dropbox/Berkeley/Project_EcoG_LPFC'
  
  # Create image directory if it doesn't exist
  if(pdf) {
    if(!dir.exists(paste(lpfc_path,'/images/stepwise_model/',freq_band,sep=''))) {
      dir.create(paste(lpfc_path,'/images/stepwise_model/',freq_band,sep=''),recursive=T)
    }
  }
  
  # Load dependencies
  source(paste(lpfc_path,'/scripts/functions/ft_single_elec_perm_pval.r',sep=''))
  source(paste(path,'/scripts/functions/single_elec_loc.r',sep=''))
  source(paste(path,'/scripts/functions/lmp.r',sep=''))
  source(paste(path,'/scripts/functions/ofc_to_subj.r',sep=''))
  source(paste(path,'/scripts/functions/get_reg_matrix.r',sep=''))
  source(paste(path,'/scripts/functions/abbrev.r',sep=''))
  rotate <- function(x) t(apply(x, 2, rev))
  print(paste('Stepwise regression analysis ',SUBID, ', elec ',elec_name, ' ', epoch,sep=''))
  
  # Avoid bad electrodes and placeholder pdf file if bad
  subj_globals = readMat(paste(path,'/data/',SUBID,'/subj_globals.mat',sep=''))
  # if(!elec_subj %in% subj_globals$good.elecs) {
  #   print(paste('Electrode',elec_subj,'is bad or discarded - skipping'))
  #   if(pdf) {
  #     print(paste('Saving pdf file /.../Project_EcoG/images/',SUBID,'/model_comp/',epoch,'_e',sprintf('%03d',elec_subj),'_',freq_band,'.pdf',sep=''))
  #     pdf(paste(path,'/images/',SUBID,'/model_comp/',epoch,'_e',sprintf('%03d',elec_subj),'_',freq_band,'.pdf',sep=''),height=10,width=7)
  #     plot(0,0,axes=F,type='n',xlab='',ylab='')
  #     text(0,0,paste(elec_subj,' is bad or discarded - skip',sep=''))
  #     dev.off()
  #   }
  #   return()
  # }
  
  # Initialize model_info
  model_info = list()
  model_info$SUBID = SUBID                          # SUBID
  model_info$elec_name = elec_name
  
  # Get anatomical location and label
  atlas = read.csv(paste(path,'/recons/all_atlas_classif.csv',sep=''))
  atlas_idx = which(atlas$SUBID %in% SUBID & atlas$elec_name %in% elec_name)
  model_info$atlas_loc = atlas$AAL[atlas_idx]
  model_info$anat_classif = atlas$anat_classif[atlas_idx] # So far limited to OFC/LPFC/Precentral/Postcentral

  # Checks for necessary files
  lpfc_path = '/Users/ignaciosaez/Dropbox/Berkeley/Project_EcoG_LPFC'
  if(!file.exists(paste(lpfc_path,'/',SUBID,'/parsed/ft_parsed_',epoch,'_',elec_name,'.mat',sep=''))) {
    print(paste('File ~/.../Project_EcoG_LPFC/',SUBID,'/parsed/ft_parsed_',epoch,'_',elec_name,'.mat does not exist',sep=''))
    print('Run ft_to_single_elec.m first!')
    return()
  }

  # Read parsed windowed data from single_elec_parse.m
  events = readMat(paste(lpfc_path, '/',SUBID,'/parsed/ft_parsed_',epoch,'_',elec_name,'.mat',sep=''))
  events = events$window.events[,,1]
  eval(parse(text=paste('window.events = events$',tolower(freq_band),sep='')))
  parse.params = as.data.frame(events$parse.params)
  nwin = dim(window.events)[2]
  
  # NEEDS TO BE FIXED
  #eval(parse(text=paste('center_window = events$parse.params$',epoch,'.center.window',sep='')))
  center_window = 21
  
  ## 1. Independent linear regressions ####
  print('Running independent linear regressions and permutation tests')
  lm_rsq  = array(data=NA,dim=c(length(reg_list),nwin))
  lm_pval = array(data=NA,dim=c(length(reg_list),nwin))
  reg_matrix = get_reg_matrix(SUBID, reg_list) # Build the regressor matrix
  if(SUBID == 's32_IR55') { reg_matrix = reg_matrix[,1:135] } # Hate to do this but it's the easiest - last trial has NaNs and is thrown out during FT preprocessing
  # Run linear regression for each regressor across all timepoints
  for(h in 1:length(reg_list)) {
    for (i in 1:nwin) {
      peak = window.events[,i]
      temp_lm = lm(peak ~ reg_matrix[h,])
      lm_rsq[h,i] = summary(temp_lm)$r.squared
      lm_pval[h,i] = -log10(summary(temp_lm)$coefficients[2,4])
    }
  }
  rownames(lm_rsq) = reg_list
  rownames(lm_pval) = reg_list
  rownames(reg_matrix) = reg_list
  
  ## 2. Select regressors that survive permutation tests and multiple comparisons ####
  # a. Subset list to regressors with any point with pval <0.05 only
  perm_pval = rep(0,length(reg_list))
  max_pval = apply(lm_pval,1,max)
  max_rsq = apply(lm_rsq,1,max)
  lm_rsq_sig = lm_rsq[which(max_pval > 1.3),]
  lm_pval_sig = lm_pval[which(max_pval > 1.3),]
  reg_matrix_sig = reg_matrix[which(max_pval > 1.3),]
  perm_pval[which(max_pval < 1.3)] = 1
  
  # b. Further subset to regressors that survive permutation test and mult comp correct
  sig_reg = rownames(lm_rsq_sig)
  if(length(sig_reg) < 1) {
    corr_pval = 1
  } else {
    corr_pval = rep(0,length(sig_reg))
    for(i in 1:length(sig_reg)) {
      reg = sig_reg[i]
      corr_pval[i] = ft_single_elec_perm_pval(SUBID, elec_name, freq_band = freq_band, reg = reg, niter = niter, epoch = epoch)
      corr_pval[i] = corr_pval[i] * corr_nelecs # correct by #elecs
    }
  }
  perm_pval[which(max_pval > 1.3)] = corr_pval
  perm_pval[perm_pval > 1] = 1 # Cap max perm pval at 1, regardless of # corrections
  
  # c. Further subset to regressors that survive multiple comparisons
  if(min(dim(as.matrix(lm_rsq_sig))) < 2) {
    lm_rsq_sig_perm = lm_rsq_sig[corr_pval < 0.05]
    lm_pval_sig_perm = lm_pval_sig[corr_pval < 0.05]
    reg_matrix_sig_perm = reg_matrix_sig[corr_pval < 0.05]
  } else {
    lm_rsq_sig_perm = lm_rsq_sig[corr_pval < 0.05,]
    lm_pval_sig_perm = lm_pval_sig[corr_pval < 0.05,]
    reg_matrix_sig_perm = reg_matrix_sig[corr_pval < 0.05,]
  }
  reg_names = rownames(lm_rsq_sig)[corr_pval < 0.05]
  nsig_reg = sum(corr_pval < 0.05)
  if(nsig_reg == 1) {
    model_comp = rownames(lm_rsq_sig)[corr_pval < 0.05] # Save name for plot
  }
  
  # d. Stop if no significant regressors
  if(nsig_reg < 1) {
    print(paste('No significant regressors for ', SUBID,' ',elec_name,sep=''))
    model_info$SUBID = SUBID
    #model_info$elec_subj = elec_subj
    model_info$ncomp = 0
    model_info$comp_names = ''
    model_info$lm_rsq_all = lm_rsq                        # Rsq for all regressors
    model_info$lm_pval_all = lm_pval                      # pval
    model_info$stepwise = numeric(0)             # Final model fit
    model_info$residual_fits = numeric(0)          # Residual fits for stepwise regressors
    model_info$reg_pval_perm = rep(1,length(reg_list))              # corrected pvals for regressor selection (for all regressors that initially survived mult.comp.) - can be >1 due to correction
    if(pdf) {
      print(paste('Saving pdf file /.../Project_EcoG_LPFC/images/',SUBID,'/stepwise_model/',epoch,'_',elec_name,'_',freq_band,'.pdf',sep=''))
      pdf(paste(lpfc_path,'/images/stepwise_model/',freq_band,'/',SUBID,'_',epoch,'_',elec_name,'_',freq_band,'.pdf',sep=''),height=10,width=7)
      plot(0,0,axes=F,type='n',xlab='',ylab='')
      text(0,0,paste('No significant regressors for ', SUBID,' ',elec_name,sep=''))
      dev.off()
    }
    # Save model results
    if(!dir.exists(paste(lpfc_path,'/',SUBID,'/model_comp/',sep=''))) {
      dir.create(paste(lpfc_path,'/',SUBID,'/model_comp/',sep=''))
    }
    save(model_info,file=paste(lpfc_path,'/',SUBID,'/stepwise_model/',freq_band,'_',epoch,'_stepwise_model_',elec_name,'.RData',sep=''))
    
    # Return results
    if(return_results) {
      return(model_info)
    } else {
      return()
    }
  }# else {
  #  print(paste(reg_subset[which(1:nsig_reg %in% reg_subset)],' survive n=',corr_nelecs,'mult comp',sep=''))
  #}
  
  ## 3. Select best regressor order ####
  # a. First, select regressor with max RSq peak
  if(dim(as.matrix(lm_rsq_sig_perm))[2] == 1) {
    rsq_max = max(lm_rsq_sig_perm)
  } else {
    rsq_max = apply(lm_rsq_sig_perm,1,max) 
  }
  reg_subset = as.numeric(which(rsq_max == max(rsq_max))) # Initialize # comps in model
  
  # b. Iterate over all possible 2nd regressors to find which one improves the fit the most; add to model
  # Keep adding until all are ordered
  if(nsig_reg > 1) { # Only if 2+ regressors are selected
    
    model_comp_fstat = model_comp_pval = matrix(NaN,nrow=nsig_reg-1,ncol=nwin) # Pairwise model comparisons # Used to include model_comp_rsq - obsolete?
    fstat_stretch = fstat_thr = matrix(NaN,nsig_reg-1)
    unchosen_regs = which(!(1:nsig_reg %in% reg_subset))
    while(length(unchosen_regs) > 0) {
      model_comp_fstat_temp = matrix(NaN,nrow=length(unchosen_regs),ncol=nwin)
      for (i in 1:nwin) {
        peak = window.events[,i]       # HG peaks at window
        
        # Simpler model
        reg_matrix_temp_1 = as.data.frame(t(rbind(peak,reg_matrix_sig_perm[reg_subset,])))  # (nreg) regression matrix: y = HG , X = (nreg)
        temp_lm_1 = lm(peak ~ ., data=reg_matrix_temp_1)                                    # Run regression and get RSq
        
        # Iterate over all possible complex models
        for (j in 1:length(unchosen_regs)) { # For all components not yet in the model
          reg_matrix_temp_2 = as.data.frame(t(rbind(peak,reg_matrix_sig_perm[c(reg_subset,unchosen_regs[j]),])))        # (nreg+1) regression matrix: y = HG , X = (nreg+1)
          temp_lm_2 = lm(peak ~ ., data=reg_matrix_temp_2)                                  # Run regression and get RSq
          anova_temp = anova(temp_lm_1,temp_lm_2)
          if(anova_temp$Df[2] == 0) { # Identical models
            model_comp_fstat_temp[j,i] = 0
          } else {
            model_comp_fstat_temp[j,i] = anova_temp$F[2]
          }
        }
      }
      # Select the best possible next regressor based on summary Fstat
      fstat_thr = qf(0.95, df1=max(anova_temp$Df[2],1),df2=anova_temp$Res.Df[2]) # Find out the f-stat threshold for p = 0.05 given the ANOVA DFs; max(df,1) avoids NaN when df=0 for comparing identical models (can happen with last few regressors)
      fstat_stretch_temp = matrix(NaN,length(unchosen_regs))
      for(h in 1:length(unchosen_regs)) {
        stretch = model_comp_fstat_temp[h,] > fstat_thr
        indices = stretch_start_end(stretch)
        if(is.na(indices[1])) {
          fstat_stretch_temp[h] = max(model_comp_fstat_temp[h,])
        } else {
          fstat_stretch_temp[h] = sum(model_comp_fstat_temp[h,indices[1]:indices[2]]) # Summary stat
        }
      }
      selected_unchosen_reg = unchosen_regs[which(fstat_stretch_temp == max(fstat_stretch_temp,na.rm=T))]
      reg_subset = c(reg_subset,selected_unchosen_reg)
      unchosen_regs = which(!(1:nsig_reg %in% reg_subset))
    }
    lm_rsq_sig_perm_ord = lm_rsq_sig_perm[reg_subset,]
    lm_pval_sig_perm_ord = lm_pval_sig_perm[reg_subset,]
    reg_matrix_sig_perm_ord = reg_matrix_sig_perm[reg_subset,]
  } else {
    lm_rsq_sig_perm_ord = lm_rsq_sig_perm
    lm_pval_sig_perm_ord = lm_pval_sig_perm
    reg_matrix_sig_perm_ord = reg_matrix_sig_perm
  }
  
  ## 4. Run incremental regressions ####
  print('Running stepwise regressions')
  # Uses a permutation strategy: for each model comparison (i.e., if 4 regressors survive, there will be 3 model comparisons):
  #   1. ANOVA Fstat timecourse is calculated for simple-complex model comparison
  #   2. Summary Fstat is calculated = sum of Fstats over fstat_thr
  #   3. A Null is calculated by shuffling labels niter times
  #   4. Calculate permuation pvals
  #   5. Repeat for each pairwise model comparison
  #
  # Full model
  # Run pairwise comparisons between (nreg) and (nreg+1) model fits
  # Store resulting F-stats and create timecourse F-stat statistic
  if(nsig_reg > 1) { # Only if 2+ regressors are selected
    model_comp_fstat = matrix(NaN,nrow=nsig_reg-1,ncol=nwin) # Pairwise model comparisons
    model_rsq = matrix(NaN,nrow=nsig_reg,ncol=nwin)
    fstat_stretch = fstat_thr = matrix(NaN,nsig_reg-1)
    
    for(h in 2:nsig_reg) {
      for (i in 1:nwin) {
        peak = window.events[,i]       # HG peaks at window
        
        # Simpler model
        reg_matrix_temp_1 = as.data.frame(t(rbind(peak,reg_matrix_sig_perm_ord[1:(h-1),])))  # (nreg) regression matrix: y = HG , X = (nreg)
        temp_lm_1 = lm(peak ~ ., data=reg_matrix_temp_1)                                     # Run regression and get RSq
        model_rsq[h-1,i]  = summary(temp_lm_1)$r.square
        
        # Complex model
        reg_matrix_temp_2 = as.data.frame(t(rbind(peak,reg_matrix_sig_perm_ord[1:h,])))        # (nreg+1) regression matrix: y = HG , X = (nreg+1)
        temp_lm_2 = lm(peak ~ ., data=reg_matrix_temp_2)                                       # Run regression and get RSq
        model_rsq[h,i]  = summary(temp_lm_2)$r.square
        
        # Compare and store F-stat timecourse
        anova_temp = anova(temp_lm_1,temp_lm_2)
        # model_comp_pval[h,i] = -log10(anova_temp$`Pr(>F)`[2])
        model_comp_fstat[h-1,i] = anova_temp$F[2]                                            # Store Fstat for model comparisons
        fstat_thr[h-1] = qf(0.95, df1=max(anova_temp$Df[2],1),df2=anova_temp$Res.Df[2]) # Find out the f-stat threshold for p = 0.05 given the ANOVA DFs; max(df,1) avoids NaN when df=0 for comparing identical models (can happen with last few regressors)
      }
    }
    # Create sum-of-F-stats statistic
    fstat_thr = mean(fstat_thr,na.rm=T) # If the last ANOVA fails (i.e. both last models are identical, which is possible), fstat_thr = NaN; solved this way (very little variation across models, <0.1%)
    for(h in 1:dim(model_comp_fstat)[1]) {
      stretch = model_comp_fstat[h,] > fstat_thr
      indices = stretch_start_end(stretch)
      if(is.na(indices[1])) {
        fstat_stretch[h] = max(model_comp_fstat[h,])
      } else {
        fstat_stretch[h] = sum(model_comp_fstat[h,indices[1]:indices[2]]) # Summary stat
      }
    }
  }
  
  # Generate Null using the same strategy but permuting labels
  if(nsig_reg > 1) { # Only if 2+ regressors are selected
    model_comp_fstat_shuf = matrix(NaN,nrow=nsig_reg-1,ncol=nwin) # Pairwise model comparisons
    fstat_stretch_perm = fstat_thr_perm = matrix(NaN,nsig_reg-1,niter)
    for(g in 1:niter) {
      for(h in 2:nsig_reg) {
        for (i in 1:nwin) {
          peak = window.events[,i]       # HG peaks at window
          
          # Simpler model
          reg_matrix_temp_1 = as.data.frame(t(rbind(peak,reg_matrix_sig_perm_ord[1:(h-1),])))  # (nreg) regression matrix: y = HG , X = (nreg)
          rnd_labels = sample(1:dim(reg_matrix_temp_1)[1])                                     # Shuffle behavioral labels (not HG!)
          reg_matrix_temp_1[,2:dim(reg_matrix_temp_1)[2]] = reg_matrix_temp_1[rnd_labels,2:dim(reg_matrix_temp_1)[2]]
          temp_lm_1 = lm(peak ~ ., data=reg_matrix_temp_1)                                     # Run regression and get RSq
          
          # Complex model
          reg_matrix_temp_2 = as.data.frame(t(rbind(peak,reg_matrix_sig_perm_ord[1:h,])))        # (nreg+1) regression matrix: y = HG , X = (nreg+1)
          rnd_labels = sample(1:dim(reg_matrix_temp_2)[1])
          reg_matrix_temp_2[,2:dim(reg_matrix_temp_2)[2]] = reg_matrix_temp_2[rnd_labels,2:dim(reg_matrix_temp_2)[2]]
          temp_lm_2 = lm(peak ~ ., data=reg_matrix_temp_2)                                       # Run regression and get RSq
          
          # Compare and store F-stat timecourse
          anova_temp = anova(temp_lm_1,temp_lm_2)
          model_comp_fstat_shuf[h-1,i] = anova_temp$F[2]                                            # Store Fstat for model comparisons
          fstat_thr_perm[h-1] = qf(0.95, df1=max(anova_temp$Df[2],1),df2=anova_temp$Res.Df[2]) # Find out the f-stat threshold for p = 0.05 given the ANOVA DFs; max(df,1) avoids NaN when df=0 for comparing identical models (can happen with last few regressors)
        }
      }
      model_comp_fstat_shuf[is.na(model_comp_fstat_shuf)] = 0 # If the complex model is a worse fit, ANOVA$F = NaN
      for(h in 1:dim(model_comp_fstat_shuf)[1]) {
        # Create sum-of-F-stats statistic
        fstat_thr_perm = mean(fstat_thr,na.rm=T) # If the last ANOVA fails (i.e. both last models are identical, which is possible), fstat_thr = NaN; solved this way (very little variation across models, <0.1%)
        stretch = model_comp_fstat_shuf[h-1,] > fstat_thr
        indices = stretch_start_end(stretch)
        if(is.na(indices[1])) {
          fstat_stretch_perm[h,g] = max(model_comp_fstat_shuf[h,])
        } else {
          fstat_stretch_perm[h,g] = sum(model_comp_fstat_shuf[h,indices[1]:indices[2]]) # Summary stat
        }
      }
    }
    
    # Permutation pvals
    pval_perm = NaN
    for (i in 1:length(fstat_stretch)) {
      pval_perm[i] = sum(fstat_stretch[i] < fstat_stretch_perm[i,])/niter
    }
    pval_perm[is.na(pval_perm)] = 1 # If two models are identical NAs are generated - fixed here
    
    ## 5. Keep the components for the most complex model that improves fit (ANOVA p<0.05 anywhere) ####
    # needs to handle no model improvements (ncomp = 1)
    # needs to stop at first no model improvement (there may be several jumps)
    if(sum(pval_perm<0.05,na.rm=T) == (nsig_reg-1)) { # Catch when all regressors are added to the model
      model_ncomp = nsig_reg
    } else {
      model_ncomp = 1
      while(pval_perm[model_ncomp] < 0.05) {
        model_ncomp = model_ncomp + 1
      }
    }
    final_model_rsq = as.matrix(lm_rsq_sig_perm_ord[1:model_ncomp,])
    final_model_pval = as.matrix(lm_pval_sig_perm_ord[1:model_ncomp,])
    if(model_ncomp == 1) {
      model_comp = rownames(lm_rsq_sig_perm_ord)[1]
    } else {
      model_comp = rownames(final_model_rsq) 
    }
    
    # Stepwise model
    stepwise_model = model_rsq[1:model_ncomp,]
    residual_fits = stepwise_model
    if(dim(final_model_rsq)[2] > 1) {
      for(i in 2:dim(stepwise_model)[1]) {
        residual_fits[i,] = stepwise_model[i,] - stepwise_model[(i-1),]
      }
    }
  }
  
  ## 5. Plot ####
  print('Plotting and saving')
  if(pdf) {
    print(paste('Saving pdf file /.../Project_EcoG_LPFC/images/stepwise_model/',freq_band,'/',SUBID,'_',epoch,'_',elec_name,'_',freq_band,'.pdf',sep=''))
    pdf(paste(lpfc_path,'/images/stepwise_model/',freq_band,'/',SUBID,'_',epoch,'_',elec_name,'_',freq_band,'.pdf',sep=''),height=10,width=7)
    col_models =  suppressWarnings(rev(brewer.pal(nsig_reg,'OrRd')))
    legend_short = abbrev(model_comp)
    
    par(mar=c(3.1,10.1,2.1,8.1),mfcol=c(5,1))
    # First plot - all RSq, independent linear regressions
    image(rotate(lm_rsq),xlab='',ylab='',yaxt='n',xaxt='n',col=heat.colors(25),main=paste(toupper(model_info$anat_classif), model_info$atlas_loc, model_info$elec_name, epoch))
    axis(side=2,at=seq(1,0,-1/(length(reg_list)-1)),labels=reg_list,las=1)
    image.plot(rotate(lm_rsq), legend.only=TRUE,col=heat.colors(25),add=T)
    
    # Second plot - RSq, significant regressors, sorted
    if(nsig_reg == 1) {
      image(as.matrix(lm_rsq_sig_perm_ord),xlab='',ylab='',yaxt='n',xaxt='n',col=heat.colors(25),main='Mult comp corrected, Fstat sorted')
      axis(side=2,at=0,labels=reg_names,las=1)
      image.plot(as.matrix(lm_rsq_sig_perm_ord), legend.only=TRUE,col=heat.colors(25),add=T) 
    } else {
      image(rotate(lm_rsq_sig_perm_ord),xlab='',ylab='',yaxt='n',xaxt='n',col=heat.colors(25),main='Mult comp corrected, Fstat sorted')
      axis(side=2,at=seq(1,0,-1/(dim(lm_rsq_sig_perm_ord)[1]-1)),labels=rownames(lm_rsq_sig_perm_ord),las=1)
      image.plot(rotate(lm_rsq_sig_perm_ord), legend.only=TRUE,col=heat.colors(25),add=T) 
    }
    
    if(nsig_reg > 1) {
      # Third plot - final number of components
      if (dim(final_model_rsq)[2] == 1) {
        image(final_model_rsq,xlab='',ylab='',yaxt='n',xaxt='n',col=heat.colors(25),main='Stepwise model regressors, Fstat sorted')
      } else {
        image(rotate(final_model_rsq),xlab='',ylab='',yaxt='n',xaxt='n',col=heat.colors(25),main='Final stepwise model regressors, Fstat sorted')
      }
      if(length(model_comp) > 1) { axis(side=2,at=seq(1,0,-1/(length(model_comp)-1)),labels=model_comp,las=1) } else { axis(side=2,at=0.5,labels=model_comp,las=1)  }
      image.plot(rotate(final_model_rsq), legend.only=TRUE,col=heat.colors(25),add=T)
      
      # Fourth plot - model RSq for incremental number of components
      plot(0,xlim=c(0,nwin),ylim=c(0,max(model_rsq,na.rm=T)),xlab='Time(s)',ylab='RSq',bty='l',xaxt='n',las=1,type='n',main='Stepwise models - RSq')
      #tick_spacing = 1000/parse.params$lag
      tick_spacing = 1000/50
      tick_loc = center_window + (c(-10:10) * tick_spacing)
      tick_loc = tick_loc[tick_loc > 0 & tick_loc <= dim(model_rsq)[2]]
      axis(side=1, las=1, at=tick_loc,labels=1:length(tick_loc) - which(tick_loc %in% center_window[[1]]))
      legend_ncomp = ''
      for(i in 1:(model_ncomp)) {
        lines(model_rsq[i,],col=col_models[i])
        legend_ncomp[i] = paste(legend_short[1:i],collapse='')
      }
      legend(x='topleft',legend=legend_ncomp,bty='n',col=col_models,lty='solid')
      
      # Fifth plot - pairwise model comparison Fstat
      plot(0,xlim=c(0,nwin),ylim=c(0,max(model_comp_fstat,na.rm=T)),xlab='Time(s)',ylab='ANOVA Fstat',bty='l',xaxt='n',las=1,type='n',main='Stepwise model comparison')
      arrows(-1,fstat_thr,100,fstat_thr,col='grey',lty='dashed')
      tick_loc = center_window + (c(-10:10) * tick_spacing)
      tick_loc = tick_loc[tick_loc > 0 & tick_loc <= dim(model_rsq)[2]]
      axis(side=1, las=1, at=tick_loc,labels=1:length(tick_loc) - which(tick_loc %in% center_window[[1]]))
      legend_ncomp = ''
      if(model_ncomp > 1) {
        for(i in 1:(model_ncomp-1)) {
          lines(model_comp_fstat[i,],col=col_models[i])
          legend_ncomp[i] = paste('(',paste(legend_short[1:i],collapse=''),') vs (',paste(legend_short[1:(i+1)],collapse=''),')',sep='')
        }
        legend(x='topleft',legend=legend_ncomp,bty='n',col=col_models,lty='solid')
      }
    }
    dev.off()
  }
  
  ## 6. Save data and return it as a list
  # Interesting data to retrieve:
  # 1. Identities of selected components: names and positions in reg_list
  # 2. Max RSq of final model
  # 3. Incremental RSq improvement, by component
  model_info$ncomp = length(model_comp)             # Number of components in final model
  model_info$comp_names = model_comp                # Names of components in final model
  model_info$comp_pos = rep(0,length(reg_list))     # Location of final components in reg_list passed
  model_info$comp_pos[match(model_info$comp_names, reg_list)] = 1:length(model_info$comp_names)
  model_info$lm_rsq = lm_rsq[match(model_info$comp_names, reg_list),]     # RSq for individual linear regressions for kept regressors only
  model_info$lm_pval = lm_pval[match(model_info$comp_names, reg_list),]   # pval
  model_info$lm_rsq_all = lm_rsq                    # Rsq for all regressors
  model_info$lm_pval_all = lm_pval                  # pval
  model_info$reg_pval_perm = perm_pval              # corrected pvals for regressor selection (for all regressors that initially survived mult.comp.) - can be >1 due to correction
  model_info$parse.params = events$parse.params
  if (nsig_reg > 1) {
    model_info$model_comp_pval_perm = pval_perm     # pvals for model comparison 
    model_info$stepwise_model = stepwise_model      # Final model fit
    model_info$residual_fits = residual_fits        # Residual fits for stepwise regressors, ordered exactly as stepwise_model
  } else {
    model_info$stepwise_model =  model_info$lm_rsq  # Final model fit
    model_info$residual_fits =  model_info$lm_rsq   # Residual fits for stepwise regressors, ordered exactly as stepwise_model
  }
  
  # Save model results
  if(!dir.exists(paste(lpfc_path,'/',SUBID,'/stepwise_model/',sep=''))) {
    dir.create(paste(lpfc_path,'/',SUBID,'/stepwise_model/',sep=''))
  }
  save(model_info,file=paste(lpfc_path,'/',SUBID,'/stepwise_model/ft_stepwise_',freq_band,'_',epoch,'_',elec_name,'.RData',sep=''))
  
  # Return results
  if(return_results) {
    return(model_info)
  }
  
}

# model_info$
# lm_rsq = results of individual linear regressions for kept regressors only            - Plot 1
# lm_rsq_sig = any p<0.05
# lm_rsq_sig_perm = survive multiple elec correction (/nelecs)
# lm_rsq_sig_perm_ord = ordered according to model improvement (based on summary Fstat) - Plot 2
# final_model_rsq = kept regressors from lm_rsq_sig_perm_ord                            - Plot 3
# model_rsq = omnibus Rsq for incremental number of regressors (nReg x nWin)            - Plot 4
# model_comp_fstat = pairwise model comparison                                          - Plot 5
