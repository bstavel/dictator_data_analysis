# Run stepwise regression for all patients with OFC-LPFC coverage
# Data is preprocessed in FT using within-grid CAR

# Prep
library(R.matlab)
source('~/Dropbox/Berkeley/Project_ECoG_LPFC/scripts/functions/ft_single_elec_stepwise.r')
reg_list = c('gamble_ind','risk','winprob','exputil','win_ind','loss_ind','rpe','regret','previous_gamble_ind','previous_risk','previous_winprob','previous_exputil','previous_win_ind','previous_loss_ind','previous_rpe','previous_regret')
SUBID_list = c('s03_AMC046','s07_GP51','s14_IR20','s25_IR38','s32_IR55')

# Run regressions for good elecs only
for(i in 1:length(SUBID_list)) {
  SUBID = SUBID_list[i]
  subj_globals = readMat(paste('~/Dropbox/Berkeley/Project_ECoG/data/',SUBID,'/subj_globals.mat',sep=''))
  good_elecs = as.numeric(subj_globals$good.elecs)
  header = readMat(paste('~/Dropbox/Berkeley/Project_ECoG/data/',SUBID,'/signal_header.mat',sep=''))
  if(i>3){
    elecs = as.character(unlist(header$header[[7]]))
  } else {
    elecs = as.character(unlist(header$header))
  }
  good_elec_names = elecs[good_elecs]
  for(j in 1:length(good_elecs)) {
    elec_name = good_elec_names[j]
    print(paste('Processing ', SUBID, ', elec ', elec_name,sep=''))
    ft_single_elec_stepwise(SUBID=SUBID,elec_name=elec_name,epoch='choice',reg_list=reg_list)
    ft_single_elec_stepwise(SUBID=SUBID,elec_name=elec_name,epoch='outcome',reg_list=reg_list)
  }
}