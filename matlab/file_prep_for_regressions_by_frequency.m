function file_prep_for_regressions_by_frequency(input_file, sub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function loads input file, reaarranges for format appropriate for R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data %
load(input_file)

% get electrode names %
elec_table = cell2table(TFwave.label);
num_elecs = size(TFwave.label, 1) ;
elec_index = 1:num_elecs ;
elec_table.index = transpose(elec_index) ;

% time indices of interest %
pre_trial_time = -.2
post_trial_time = 2

% get num of trials %
nTrials = size(TFwave.trialinfo(:, 1), 1);

%% theta %%
low_freq = 4;
high_freq = 8;

% convert to specific grequency %
freq_index = (TFwave.freq >= low_freq & TFwave.freq < high_freq);
freq_df = TFwave.powspctrm(1:nTrials, :, freq_index, :);
freq_average = mean(freq_df, 3);
theta_average_squeeze = squeeze(freq_average);

% concactenate and cut hg %
for idx = 1:nTrials
  tIdx = TFwave.trialinfo(idx, 1); % since some trials are skipped, need to reindex
  % get indices between second 0 and 3 where, 0 is presentation time %
  indices_of_interest = (TFwave.time < post_trial_time & TFwave.time > pre_trial_time) ; % because there was iti included, and -.2 + .6 = .4
  % cut the extrad padding on each trial window %
  temp_hg =  squeeze(theta_average_squeeze(idx, 1:num_elecs, indices_of_interest)) ;
  % sanity check to save elecs order %
  elecs_counter = size(temp_hg, 2);
  temp_hg(:, (elecs_counter + 1)) = 1:num_elecs ;
  temp_hg(:, (elecs_counter + 2)) = tIdx ;
  % concactenate across trials $
  if idx == 1
    hg_prepped = temp_hg ;
  else
    hg_prepped = vertcat(hg_prepped, temp_hg) ;
  end

end


% save data %
csvwrite(sprintf('~/Projects/dictator_analysis/dictator_game/dg_behave_analysis/munge/%s_low-%d_high-%d_munge_presentation_locked_extended.csv', sub, low_freq, high_freq), hg_prepped)

%% alpha %%
low_freq = 8;
high_freq = 12;

% convert to specific grequency %
freq_index = (TFwave.freq >= low_freq & TFwave.freq < high_freq);
freq_df = TFwave.powspctrm(1:nTrials, :, freq_index, :);
freq_average = mean(freq_df, 3);
theta_average_squeeze = squeeze(freq_average);

% concactenate and cut hg %
for idx = 1:nTrials
  tIdx = TFwave.trialinfo(idx, 1); % since some trials are skipped, need to reindex
  % get indices between second 0 and 3 where, 0 is presentation time %
  indices_of_interest = (TFwave.time < post_trial_time & TFwave.time > pre_trial_time) ; % because there was iti included, and -.2 + .6 = .4
  % cut the extrad padding on each trial window %
  temp_hg =  squeeze(theta_average_squeeze(idx, 1:num_elecs, indices_of_interest)) ;
  % sanity check to save elecs order %
  elecs_counter = size(temp_hg, 2);
  temp_hg(:, (elecs_counter + 1)) = 1:num_elecs ;
  temp_hg(:, (elecs_counter + 2)) = tIdx ;
  % concactenate across trials $
  if idx == 1
    hg_prepped = temp_hg ;
  else
    hg_prepped = vertcat(hg_prepped, temp_hg) ;
  end

end


% save data %
csvwrite(sprintf('~/Projects/dictator_analysis/dictator_game/dg_behave_analysis/munge/%s_low-%d_high-%d_munge_presentation_locked_extended.csv', sub, low_freq, high_freq), hg_prepped)


%% beta %%
low_freq = 12;
high_freq = 30;

% convert to specific grequency %
freq_index = (TFwave.freq >= low_freq & TFwave.freq < high_freq);
freq_df = TFwave.powspctrm(1:nTrials, :, freq_index, :);
freq_average = mean(freq_df, 3);
theta_average_squeeze = squeeze(freq_average);

% concactenate and cut hg %
for idx = 1:nTrials
  tIdx = TFwave.trialinfo(idx, 1); % since some trials are skipped, need to reindex
  % get indices between second 0 and 3 where, 0 is presentation time %
  indices_of_interest = (TFwave.time < post_trial_time & TFwave.time > pre_trial_time) ; % because there was iti included, and -.2 + .6 = .4
  % cut the extrad padding on each trial window %
  temp_hg =  squeeze(theta_average_squeeze(idx, 1:num_elecs, indices_of_interest)) ;
  % sanity check to save elecs order %
  elecs_counter = size(temp_hg, 2);
  temp_hg(:, (elecs_counter + 1)) = 1:num_elecs ;
  temp_hg(:, (elecs_counter + 2)) = tIdx ;
  % concactenate across trials $
  if idx == 1
    hg_prepped = temp_hg ;
  else
    hg_prepped = vertcat(hg_prepped, temp_hg) ;
  end

end


% save data %
csvwrite(sprintf('~/Projects/dictator_analysis/dictator_game/dg_behave_analysis/munge/%s_low-%d_high-%d_munge_presentation_locked_extended.csv', sub, low_freq, high_freq), hg_prepped)

%% HFA %%
low_freq = 70;
high_freq = 200;

% convert to specific grequency %
freq_index = (TFwave.freq >= low_freq & TFwave.freq < high_freq);
freq_df = TFwave.powspctrm(1:nTrials, :, freq_index, :);
freq_average = mean(freq_df, 3);
theta_average_squeeze = squeeze(freq_average);

% concactenate and cut hg %
for idx = 1:nTrials
  tIdx = TFwave.trialinfo(idx, 1); % since some trials are skipped, need to reindex
  % get indices between second 0 and 3 where, 0 is presentation time %
  indices_of_interest = (TFwave.time < post_trial_time & TFwave.time > pre_trial_time) ; % because there was iti included, and -.2 + .6 = .4
  % cut the extrad padding on each trial window %
  temp_hg =  squeeze(theta_average_squeeze(idx, 1:num_elecs, indices_of_interest)) ;
  % sanity check to save elecs order %
  elecs_counter = size(temp_hg, 2);
  temp_hg(:, (elecs_counter + 1)) = 1:num_elecs ;
  temp_hg(:, (elecs_counter + 2)) = tIdx ;
  % concactenate across trials $
  if idx == 1
    hg_prepped = temp_hg ;
  else
    hg_prepped = vertcat(hg_prepped, temp_hg) ;
  end

end


% save data %
csvwrite(sprintf('~/Projects/dictator_analysis/dictator_game/dg_behave_analysis/munge/%s_low-%d_high-%d_munge_presentation_locked_extended.csv', sub, low_freq, high_freq), hg_prepped)


writetable(elec_table, sprintf('~/Projects/dictator_analysis/dictator_game/dg_behave_analysis/munge/%s_electrodes_presentation_TFwaves.csv', sub))

return
