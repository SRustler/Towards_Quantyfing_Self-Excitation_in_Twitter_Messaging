%================================
%=  Rustler Stefan, 2014        =
%=  <rustlers@student.ethz.ch>  =
%================================
%
% Data snippet:
%   presC{1}(i) - MT-id|presC{2}(i) - tweet-id|presC{3}(i) - timestamp|presC{4}(i) - user-id|presC{5}(i) - count
%   0                  |334979275864297472    |2013-05-16 10:30:57    |54491252             |1
%   0                  |334979276325654529    |2013-05-16 10:30:57    |376712526            |0
%   334979275864297472 |334988950307487745    |2013-05-16 11:09:24    |400642052            |-1    
%
% Input:
% - presC: Cell of present data that contains MTid (if applicable), tweet-id, time stamp, user-id, #RTs (-1 if RT) 
% - tevnts: Transformed time stamps in hours. Note that numel(tevnts)=length(presC). If you want to handle original time stamps, enter 24*datenum(presC{1,3},'yyyy-mm-dd HH:MM:SS') as tevnts.
%
% Output: 
% - RTMTtau_pres: Matrix with each row containing RT-index, respective MT-index, and tau between MT and RT (in units of days).
% 
% Note: 
% - RTMTtau_pres(:,3) will be in units of tevnts, i.e. hours.

function [RTMTtau_pres] = get_taus_pres(presC,tevnts) 
    SetSize = length(presC{1});         %Size of the whole set
    % Search data for RTs
    RTMTtau_pres=zeros(0.4*SetSize,3);  %Reserve space for matrix that will contain the index of the RT (1st col), index of the MT (2nd col) and time to MT (3rd col)
    ii = 0;                             %Index for RTMTtau_pres (ii will be less than i!)
    % Go through all tweets to search for a RT marked by count = -1
    for i=1:SetSize 
        if presC{5}(i)==-1                      %If the current tweet is a RT, i.e. count (of RTs) = "-1".
            time_RT = tevnts(i);                %%ONLY DIFFERENCE TO get_taus_pres_old.m!             
            for j=i-1:-1:1                      %Scan whole history of tweets previous to RT, starting from most recent in history
                if (presC{2}(j) == presC{1}(i)) %If ID of index j equals MT-ID of RT of index i.
                    time_MT = tevnts(j);        %%ONLY DIFFERENCE TO get_taus_pres_old.m!             
                    ii = ii + 1;                %RT-counter
                    RTMTtau_pres(ii,1) = i;     %Assign index of RT w.r.t. originial set (presC)
                    RTMTtau_pres(ii,2) = j;     %Assign index of MT w.r.t. originial set (presC)
                    RTMTtau_pres(ii,3) = 1/24* (time_RT - time_MT); %This is converted to units of days!
                    break
                end 
            end
        end   
        getting_taus = [num2str((i/SetSize)*100,'%.2f'),' %'] %Check progress... percentage = ii/lenght(...): presC{1} has non-zero entries when it's a RT!

%         getting_RTs = [num2str((ii+1)/length(nonzeros(presC{1}))*100,'%.2f'),' %'] %Check progress... percentage = ii/lenght(...): presC{1} has non-zero entries when it's a RT!
    end 
    RTMTtau_pres = RTMTtau_pres(1:ii,1:3);      %Discard unoverwritten zeros such that they are not plotted in histogram!
end


