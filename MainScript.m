%================================
%=  Rustler Stefan, 2014        =
%=  <rustlers@student.ethz.ch>  =
%================================
%
%% ===== MAIN SCRIPT =====
%
% This script creates and visualizes different response time distributions. 
% For the PRESENT data it can detrend time series ('binsize','unit',['smoothingParam']), regard users of certain activity only ('num'), systematically disregard extreme response times per for each user ('cutoff'). 
% For the PAST data it can detrend time series ('binsize', [unit=!=day, 'smoothingParam']), set thresholds (string distance) for RT-MT linking ('levenshtein')
% 
% Input:
% - Present data: unibrennt_daten_present.txt, unibrennt_daten_present_onlyMTs.txt
% - Past data: unibrennt_daten_past.txt, unibrennt_daten_past_date-string-MTid.txt
% 
% Notes:
% - get_user_taus() takes especially long to run. load('/Results/MAIN_PRESENT.mat') to get the data. You can then plot different CCDFs.
% - The variables in paranthesis() are those relevant to performing the respective task in the function. Variables within brackets[] are not function arguments but would have to be altered in the function itself .
% 
clear;
addpath(genpath(pwd));  %Adds all directories and subdirectories of the present working directory (Crucial for functions calling and data reading).
%
%% ===== Present data =====

% == Declarations and file reading ==

fileID1 = fopen('unibrennt_daten_present.txt');
fileID2 = fopen('unibrennt_daten_present_onlyMTs.txt');
presdata_All = textscan(fileID1, '%f %f %s %f %d','Delimiter','\t'); %File format: MT-id, tweet-id, time stamp, user-id, #RTs
presdata_MTs = textscan(fileID2, '%f %f %s %f %d','Delimiter','\t');
presdata_MTs = presdata_MTs{3}; %Only the time stamps are needed.
fclose('all'); clear fileID1 fileID2;

binsize = 240;      %Binsize in minutes used for getting rate R(t) which is then spline-fitted. Choose rather high value for present-data due to rather low tweet activity.
unit = 'month';     %Consider longer range due to rather low tweet activity in present-data.
idx = 1:31;         %Index for ii-loop, i.e. through which units, e.g. months, you want to loop. Note that the first and last month are incomplete, nevertheless they should be included s.t. t* and thus tau* is complete.
visual = 0;         %Choose 1 if R(t), N(t) etc. shall be plotted.

activity = 'RTs';   %'tweets'; %'RTs'; %State whether an "active" user is defined via number of RTs or tweets.
num = 300;          %Minimal activity of user to be plotted in figure 2.
cutoff = 3;         %Number of last, i.e. longest, RTs (most extreme events) to be cut off when plotting and before aggregating the taus again. 


% == Calculations and plotting ==

% Detrending:
[tstars_pres,gofcell] = detrend_anyFit(binsize,unit,idx,presdata_MTs,presdata_All{3},visual); 

% Process transformed times:
[RTMTtaustar_pres] = get_taus_pres(presdata_All,tstars_pres);   %Obtain transformed response times in unit of days.
RTMTtaustar_pres(:,3) = 1/24 * RTMTtaustar_pres(:,3); %Unit: d->h
taustars_pres = RTMTtaustar_pres(:,3); %Unit: h
    fig1 = figure;
        plot_RTs(RTMTtaustar_pres,length(presdata_All{1}),'loglog')  %Four subplots visualizing RTMTtau: idxMT vs tau, idxMT vs idxRT, tau-histogram, tau-CDF.
        mtit(fig1,'Detrended RT-statistics of present');
    %fig1a (created in plot_PLwithERR)
        positaus = sort(taustars_pres); positaus = positaus(positaus>0); %Create vector of positive taus s.t. the following function works. ISSUE: There are still negative taus in the data (~10%)!
        plot_PLwithERR(positaus) %Plot CCDF with Aaron Clauset's PL codes
    %fig2 (created in plot_user_taus)
        [user_taustars] = get_user_taus(presdata_All,RTMTtaustar_pres);  %With presdata_All still containing the userIDs and RTMTtau all the taus, this function gets tau-vectors for each of the 2400 users in the present-data. Row: userID, #RTs, #tweets, taus.
        [user_taustars_AggCut] = plot_user_taus(user_taustars,activity,num,cutoff,'loglog'); %Note: user_tausstars_"Aggregated-Cutoff-Transformed" is only of those taus whose users fulfill the activity criterion!
        title(['Detrended user statistics of present(#',activity,'>',num2str(num),', last ',num2str(cutoff),' RTs per user omitted)']);

% Process original times:
[RTMTtau] = get_taus_pres(presdata_All,24*datenum(presdata_All{1,3},'yyyy-mm-dd HH:MM:SS'));  %Creates matrix with a row (idxRT, idxMT, tau) for each RT found. The number printed goes until the 19266 (total number of tweets in the present-data). Note that the second argument is contained in the first in fact. This is only to enable putting the transformed time as 2nd arg instead in other places of MainScript.m
taus = RTMTtau(:,3); %Unit: h
    fig3 = figure;
        plot_RTs(RTMTtau,length(presdata_All{1}),'loglog')  %Four subplots visualizing RTMTtau: idxMT vs tau, idxMT vs idxRT, tau-histogram, tau-CDF.
        mtit(fig3,'Undetrended RT-statistics of present');  
    %fig3a (created in plot_PLwithERR)
        positaus = sort(taus); positaus = positaus(positaus>0); %Create vector of positive taus s.t. the following function works. ISSUE: There are still negative taus in the data (~10%)!    
        plot_PLwithERR(positaus) %Plot CCDF with Aaron Clauset's PL codes
    fig4 = figure;    
        [user_taus] = get_user_taus(presdata_All,RTMTtau);  %With presdata_All still containing the userIDs and RTMTtau all the taus, this function gets tau-vectors for each of the 2400 (=USERCOUNT) users in the present-data. Row: userID, #RTs, #tweets, taus.
        [user_taus_AggCut] = plot_user_taus(user_taus,activity,num,cutoff,'loglog'); %Note: AggregUserTausWithCutoff is only of those taus whose users fulfill the activity criterion!
        title(['Undetrended user statistics of present (#',activity,'>',num2str(num),', last ',num2str(cutoff),' RTs per user omitted)']);

% Plot and compare response times:
    fig5 = figure;                          %Colors are chosen such that pure colors -> original time
        axes = 'loglog';                    %Respective slightly mixed color -> respective transformed time
        plot_taus(1/24*taus,axes,[1 0 0]);  %For plot_taus divide by 24h/d since plot_taus has axis labels in days and taus is in hours.
        hold on;
        plot_taus(1/24*taustars_pres,axes,[1 0.5 0]);  
        plot_taus(1/24*user_taus_AggCut,axes,[0 0 1]);
        plot_taus(1/24*user_taustars_AggCut,axes,[0 0.5 1]);
        legend('\tau - Original response times','\tau* - Transformed response times','\tau_{users} - Aggregated user response times','\tau*_{users} - Aggregated transformed user response times','Location','SouthWest');
        mtit(fig5,'Response time statistics of present - Comparison');
        hold off;
        
    %fig6 (created in plot_user_taus)
        %[taustars] = tstar2taustar_present(tstars,RTMTtau,'loglog'); 
        %This command gives exactly the same taustars as calculated differently in line 34. The codes are self-consistent! The command also plots the same plot as lines 60 and 62. 

save([pwd,'/Results/MAIN_PRESENT.mat']); %clear;
saveas(fig1,[pwd,'/Results/Detrended_RT-statistics_of_present.eps']); 
%saveas(fig2,'Detrended_user_statistics_of_present.eps'); happens in function 'plot_user_taus.m'
saveas(fig3,[pwd,'/Results/Undetrended_RT-statistics_of_present.eps']); 
saveas(fig4,[pwd,'/Results/Unetrended_user_statistics_of_present.eps']); 
saveas(fig5,[pwd,'/Results/Response_time_statistics_of_present-Comparison.eps']); 

%% ===== Past data =====

% == Declarations and file reading ==

fileID1 = fopen('RawData/unibrennt_daten_past.txt');
fileID2 = fopen('/RawData/unibrennt_daten_past_date-string-MTid.txt');
pastdata_All = textscan(fileID1, '%s %s %s','Delimiter','\t'); %File format: time stamp, tweet message, username
pastdata_MTs = textscan(fileID2, '%s %s %s','Delimiter','\t');
pastdata_MTs = pastdata_MTs{1}; %Only the time stamps are needed.
fclose('all'); clear fileID1 fileID2;

binsize = 10;   %Binsize in minutes used for getting rate R(t) which is then spline-fitted. Choose rather high value for present-data due to rather low tweet activity.
unit = 'day';   %Consider longer range due to rather low tweet activity in present-data.
idx = 1:33;     %Index for ii-loop, i.e. through which units, e.g. months, you want to loop. Note that the first and last date are incomplete, nevertheless they should be included s.t. t* and thus tau* is complete.
visual = 0;     %Choose 1 if R(t), N(t) etc. should be plotted for each day(!).
LevenshteinDist = 0:20:140; % Analyze with different Levenshtein distances. 0:20:140 definitely works. Errors might occur in get_taus_past for different choices.
levidx = 3;     %Manually chosen "best" Levenshtein distance for a first glance. levidx !<=! numel(LevenshteinDist). Note: This variable will be overwritten after first use.

% == Calculations and plotting ==

[RTMTtau_pastvsLeven] = get_taus_past(pastdata_All,LevenshteinDist);        %Currently commented out in function: save(['..Results/RTs_vs_Levenshtein/Data_LevenshteinDist',num2str(levenshtein),'.mat'],'RTMTtau_past')

% Process original times:
    fig_past1 = figure;
        plot_RTs( {levidx},length(pastdata_All{1}),'loglog')  %Four subplots visualizing RTMTtau: idxMT vs tau, idxMT vs idxRT, tau-histogram, tau-CCDF.        
        fig1tit = ['Undetrended RT-statistics of past for Levenshtein distance ',num2str((levidx-1)*(max(LevenshteinDist)/(numel(LevenshteinDist)-1)))];
        mtit(fig_past1,fig1tit); %Store the title for data saving at the very bottom.  
    fig_past1a = figure;
        positaus = sort(RTMTtau_pastvsLeven{levidx}(:,3)); positaus = positaus(positaus>0); %Create vector of positive taus s.t. the following function works. ISSUE: There are still negative taus in the data (~10%)!
        plot_PLwithERR(positaus) %Plot CCDF with Aaron Clauset's PL codes 
    fig_past2 = figure; 
        legendcell = cell(1,2);                                 %Preallocate some memory for cell of figure legend.
        levidx = 1;
        for leven = LevenshteinDist
            taus_past = RTMTtau_pastvsLeven{levidx}(:,3) / 24;  %Change to units of day for plot_tau.m!
%             load([pwd,'/Results/RTs_vs_Levenshtein/Data_LevenshteinDist',num2str((leven-1)*(max(LevenshteinDist)/(numel(LevenshteinDist)))),'.mat'],'RTMTtau_past'); %This line is just for testing purposes. Load the Levenshtein data instead of lengthily calculating it again.
            plot_taus(taus_past,'loglog',[(1-leven/140) leven/140 0]);  %Two subplots visualizing taus: tau-histogram, tau-CCDF (bottom two subplots of plot_RTs)
            hold on;
            legendcell{levidx} = ['Levenshtein distance = ',num2str(leven)];
            levidx = levidx + 1;
        end
        legend(legendcell,'Location','SouthWest');
        mtit(fig_past2,'Undetrended RT-statistics of past for different Levenshtein distances');   
        hold off;

% Detrending:  
[tstars_past,gofcell_past] = detrend_anyFit(binsize,unit,idx,pastdata_MTs,pastdata_All{1},visual); 

% Process transformed times: 
    fig_past3 = figure;
        legendcell = cell(1,2); %Preallocate some memory for cell of figure legend.
        levidx = 1;
        for leven = LevenshteinDist
            if leven == 0; continue; end %ISSUE: The scale of this particular one seems to be off!
            [taustars_past] = ts2taus(RTMTtau_pastvsLeven{levidx},tstars_past); %Note that this function needs RTMTtau which it loads in the function itself.
            plot_taus(taustars_past / 24,'loglog',[(1-leven/140) leven/140 0]);  %Two subplots visualizing taus: tau-histogram, tau-CCDF (bottom two subplots of plot_RTs            hold on;
            hold on;
            legendcell{1,levidx} = ['Levenshtein distance = ',num2str(leven)];
            levidx = levidx + 1;
        end
        legend(legendcell,'Location','SouthWest');
        mtit(fig_past3,'Detrended RT-statistics of past for different Levenshtein distances');
        hold off;
                                                                                                                                          
%fig_past3i (created in plot_PLwithERR for each Levenshtein distance)
        legendcell = cell(1,2); %Preallocate some memory for cell of figure legend.
        levidx = 1;
        for leven = LevenshteinDist
            if leven == 0; continue; end %ISSUE: The scale of this particular one seems to be off!
            [taustars_past] = ts2taus(RTMTtau_pastvsLeven{levidx},tstars_past); %Note that this function needs RTMTtau which it loads in the function itself.
            positaus = sort(taus_past); positaus = positaus(positaus>0); %Create vector of positive taus s.t. the following function works. ISSUE: There are still negative taus in the data (~10%)!
            plot_PLwithERR(positaus) %Plot CCDF with Aaron Clauset's PL codes
        end
        
% Plot and compare response times:
   fig_past4 = figure;        
        leven = LevenshteinDist(1);                     %First Levenshtein distance chosen
        taus_past = RTMTtau_pastvsLeven{1}(:,3) / 24;   %First Levenshtein distance chosen
        [taustars_past] = ts2taus(RTMTtau_pastvsLeven{1},tstars_past); %Note that this function needs RTMTtau which it loads in the function itself.
        plot_taus(taustars_past / 24,'loglog',[(1-leven/140) leven/140 0]);  %For plot_taus divide by 24h/d since plot_taus has axis labels in days and taus is in hours.
        hold on;
        plot_taus(taus_past / 24,'loglog',[(1-leven/140) leven/140 0]);  %For plot_taus divide by 24h/d since plot_taus has axis labels in days and taus is in hours.            
        
        leven = LevenshteinDist(end);                    %Last Levenshtein distance chosen
        taus_past = RTMTtau_pastvsLeven{end}(:,3) / 24;  %Last Levenshtein distance chosen
        [taustars_past] = ts2taus(RTMTtau_pastvsLeven{end},tstars_past); %Note that this function needs RTMTtau which it loads in the function itself.
        plot_taus(taustars_past / 24,'loglog',[leven/140 (1-leven/140) 0]);  %For plot_taus divide by 24h/d since plot_taus has axis labels in days and taus is in hours.
        hold on;
        plot_taus(taus_past / 24,'loglog',[leven/140 (1-leven/140) 0]);  %For plot_taus divide by 24h/d since plot_taus has axis labels in days and taus is in hours.            
        
        legend(['\tau - Original response times for l = ',num2str(LevenshteinDist(1))],['\tau - Original response times for l = ',num2str(LevenshteinDist(end))],['\tau* - Transformed response times for l = ',num2str(LevenshteinDist(1))],['\tau* - Transformed response times for l = ',num2str(LevenshteinDist(end))],'Location','SouthWest');
        mtit(fig_past4,'Response time statistics of present - Comparison');
        hold off;        
        
save([pwd,'/Results/MAIN_PRESENT.mat']); %clear;
saveas(fig_past1,[pwd,'/Results/',fig1tit,'.pdf']);   
saveas(fig_past2,[pwd,'/Results/Undetrended_RT-statistics_of_past_for_different_Levenshtein_distances.pdf']);
saveas(fig_past3,[pwd,'/Results/Detrended_RT-statistics_of_past_for_different_Levenshtein_distances.pdf']);
% Note: Clauset's plots, fig_past3i, are not automatically saved!
saveas(fig_past4,[pwd,'/Results/Response_time_statistics_of_past-Comparison.pdf']); 