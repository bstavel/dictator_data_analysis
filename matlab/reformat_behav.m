function reformat_behav( DATA_behav, raw_behav, photod_offset, reformat )
%make raw behavioral task data clean (1 row = 1 round), add regressors
%TO ADD MORE ROUND INFO:
%   - add col name to "headers"
%   - add calc info in for loop
%   - change .csv line writing info

% import behav data
load([DATA_behav, raw_behav])
nrows = size(data,1);

% adjust time by photo diode offset
first_choice = data{[data{:,1}]'==1 & strcmp(data(:,2),'ChoicePresentation'),4};
diff = first_choice-photod_offset;

data(:,4) = mat2cell([data{:,4}]'-diff,ones(nrows,1),1);

% for each round, extract response data
Nrounds = max(cell2mat(data(:,1)));
headers = {'round','newround.time','choice.time','buttonpress.time',...
    'reveal.time','self.payoff','other.payoff','self.foregone',...
    'other.foregone','ineq','side.chosen','RT'};
round_info = cell(0,size(headers,2));


for i=1:Nrounds
    round_info(end+1,1) = {i}; %init new round_info row
   
    R = data([data{:,1}]'==i,:); %select data rows for this round

    round_info(end,2) = R(strcmp(R(:,2),'NewRound'),4); %new round time
    round_info(end,3) = R(strcmp(R(:,2),'ChoicePresentation'),4); %choice onset time
    
    round_info(end,5) = R(strcmp(R(:,2),'Feedback'),4); %outcome reveal time
    
    response = R(strcmp(R(:,2),'SubjectConfirming'),:); %pull response row
    if size(response,1)>0 %collected response!
        button = response{1,3};
        round_info(end,11) = {button}; %which button pressed
        round_info(end,4) = response(1,4); %button press time
        
        if strcmp(button,'L') 
            round_info(end,6:9) = response(1,5:8); %save amounts as L = payoffs & R = foregone
        else
            round_info(end,6:9) = response(1,[7:8, 5:6]); %save amounts as R = payoffs & L = foregone
        end
        
        round_info(end,10) = {round_info{end,6} - round_info{end,7}}; % save inequality as self-other 
        round_info(end,12) = {round_info{end,4} - round_info{end,3}}; % RT
    else %no rseponse!
        round_info(end, [6 7 10 12]) = {[0] [0] [0] [0]};
    end
        
end
     
save([DATA_behav, reformat, '.mat'],'round_info','headers');

fid = fopen([DATA_behav, reformat, '.csv'], 'w');
fprintf(fid, [repmat(['%s,'],1,size(headers,2)-1),'%s\n'], headers{1,:});
for i=1:size(round_info,1)
    fprintf(fid, [repmat(['%f,'],1,11),'%s\n'], round_info{i,:});
end
fclose(fid);



end

