 function hg_gamma_file_prep_for_regressions(input_file, sub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function loads input file, reaarranges for format appropriate for R, should check to exclude timeouts- will check in R
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load data %
hg_raw = load(input_file)

% get electrode names and index them in order for later merging %
elec_table = cell2table(hg_raw.dataAvg2.label);
num_elecs = size(elec_table, 1) ;
elec_index = 1:num_elecs ;
elec_table.index = transpose(elec_index) ;

% get num of trials %
nTrials =  size(hg_raw.dataAvg2.trialinfo, 1)

% concactenate and cut hg %
for idx = 1:nTrials

   % get indices between second -.2 and 1 where 0 is choice time %
   pre_presentation_time = -hg_raw.dataAvg2.trialinfo(idx, 12) -.2 ;
   indices_of_interest = find(hg_raw.dataAvg2.time{idx} < 1 & hg_raw.dataAvg2.time{idx} > pre_presentation_time) ;
   % cut the extrad padding on each trial window %
   temp_hg =  hg_raw.dataAvg2.trial{idx}(1:num_elecs, indices_of_interest) ;
   % sanity check to save elecs order %
   temp_hg = [transpose(1:num_elecs) transpose(repelem(idx, num_elecs)) temp_hg];
   % concactenate across trials $
   if idx == 1
     hg_prepped = temp_hg ;
   elseif size(temp_hg, 2) > size(hg_prepped, 2)
     na_temp = nan(size(hg_prepped, 1), size(temp_hg, 2));
     tIdx1 = transpose(1:size(hg_prepped, 1)) ;
     tIdx2 = transpose(1:size(hg_prepped, 2)) ;
     na_temp(tIdx1, tIdx2) = hg_prepped;
     hg_prepped = vertcat(na_temp, temp_hg) ;
   elseif size(hg_prepped, 2) > size(temp_hg, 2)
     na_temp = nan(num_elecs, size(hg_prepped, 2));
     tIdx1 = transpose(1:size(temp_hg, 1)) ;
     tIdx2 = transpose(1:size(temp_hg, 2)) ;
     na_temp(tIdx1, tIdx2) = temp_hg;
     hg_prepped = vertcat(hg_prepped, na_temp) ;
   else
     hg_prepped = vertcat(hg_prepped, temp_hg) ;
   end

end


% save data %
csvwrite(sprintf('~/Projects/dictator_analysis/dictator_game/dg_behave_analysis/munge/%s_hg_munge_choice_locked.csv', sub), hg_prepped)
writetable(elec_table, sprintf('~/Projects/dictator_analysis/dictator_game/dg_behave_analysis/munge/%s_electrodes_choice_locked.csv', sub))

return
