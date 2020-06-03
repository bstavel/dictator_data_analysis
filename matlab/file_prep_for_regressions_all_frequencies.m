function file_prep_for_regressions_all_frequencies(input_file, sub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function loads input file, reaarranges for format appropriate for R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data %
load(input_file)

% time indices of interest %
pre_trial_time = -.2
post_trial_time = 2

% get num of trials %
nTrials = size(TFwave.trialinfo(:, 1), 1);

% get num of electrode %
elec_table = cell2table(TFwave.label);
num_elecs = size(TFwave.label, 1) ;

% pull out the data %
freq_df = TFwave.powspctrm;
freq_spec = TFwave.freq;

for fIdx  = 96:size(freq_spec, 2)
% concactenate and cut hg %
  for idx = 123:nTrials
    tIdx = TFwave.trialinfo(idx, 1); % since some trials are skipped, need to reindex
    % get indices between second 0 and 3 where, 0 is presentation time %
    indices_of_interest = (TFwave.time < post_trial_time & TFwave.time > pre_trial_time) ; % because there was iti included, and -.2 + .6 = .4
    % cut the extrad padding on each trial window %
    temp_hg =  squeeze(freq_df(idx, 1:num_elecs, fIdx, indices_of_interest)) ;
    % sanity check to save elecs order %
    elecs_counter = size(temp_hg, 2);
    temp_hg(:, (elecs_counter + 1)) = 1:num_elecs ;
    temp_hg(:, (elecs_counter + 2)) = tIdx ;
    temp_hg(:, (elecs_counter + 3)) = freq_spec(fIdx) ;
    % concactenate across trials $
    if idx == 1 && fIdx == 1
      hg_prepped = temp_hg ;
    else
      hg_prepped = vertcat(hg_prepped, temp_hg) ;
    end

  end

end


% save data %
csvwrite(sprintf('~/Projects/dictator_analysis/dictator_game/dg_behave_analysis/munge/%s_all_frequencies_munge_presentation_locked_extended.csv', sub), hg_prepped)

return
