 function hg_gamma_file_prep_for_regressions(input_file, sub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function loads input file, reaarranges for format appropriate for R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data %
hg_raw = load(input_file)

% get electrode names %
elec_table = cell2table(hg_raw.dataAvg2.label);
num_elecs = size(elec_table, 1) ;
elec_index = 1:num_elecs ;
elec_table.index = transpose(elec_index) ;

% get num of trials %
nTrials =  size(hg_raw.dataAvg2.trialinfo, 1)

% concactenate and cut hg %
for idx = 1:nTrials
   % get indices between second 0 and 3 where, 0 is presentation time %
   indices_of_interest = find(hg_raw.dataAvg2.time{idx} < 1 & hg_raw.dataAvg2.time{idx} > -.75) ;
   % cut the extrad padding on each trial window %
   temp_hg =  hg_raw.dataAvg2.trial{idx}(1:num_elecs, indices_of_interest) ;
   % sanity check to save elecs order %
   elecs_counter = size(temp_hg, 2);
   temp_hg(:, (elecs_counter + 1)) = 1:num_elecs ;
   temp_hg(:, (elecs_counter + 2)) = idx ;
   % concactenate across trials $
   if idx == 1
     hg_prepped = temp_hg ;
   else
     hg_prepped = vertcat(hg_prepped, temp_hg) ;
   end

end


% save data %
csvwrite(sprintf('~/Projects/dictator_analysis/dictator_game/dg_behave_analysis/munge/%s_hg_munge_choice_locked_cut.csv', sub), hg_prepped)
writetable(elec_table, sprintf('~/Projects/dictator_analysis/dictator_game/dg_behave_analysis/munge/%s_electrodes_choice_locked_cut.csv', sub))

return
