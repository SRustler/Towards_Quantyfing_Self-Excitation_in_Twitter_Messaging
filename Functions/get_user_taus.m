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
% - RTMTtau: Matrix with each row containing RT-index, respective MT-index, and tau between MT and RT.
%
% Output:
% - user_taus: Cell in which each row contains corresponds to one user with unique user-ID, number of active RTs, number of tweets, vector of individual responds times.
%
% Notes:
% - You can also load('TauDistr_IndividualUsers.mat','C','RTMTtau','SetSize') and only plot the results (line 40)
% - For transformed time, calculate RTMTtau_trafod first with RTs_present_trafod.m, then load('TauDistr_IndividualUsers_detrended.mat','C','SetSize')
%

function [user_taus] = get_user_taus(presC,RTMTtau)
    user = transpose(unique(presC{1,4}));   %Vector of unique user-ids
    user_taus = cell(4,length(user));       %Define cell array that will contain user-id, his/her taus, numTs and numRT
    ii = 0;
    for u = user                        %Loop through all users
        tau = zeros(1,10);              %Initially assume a maximum of 10 responses of any user (just for initialization)
        numRT = 0;                      %Count number of active RTs, i.e. the ones that the current user is tweeting (as opposed to passive RTs, i.e. tweets that are RTed by others).
        numTs = 0;                      %Count number of tweets of current user.
        for i = 1:length(presC{1})      %For each current user scan whole set.
            if (presC{4}(i)==u)         %If some tweet of current user is found...
                numTs = numTs + 1;      
                if (presC{5}(i)==-1)    %If tweet is a RT...
                    numRT = numRT + 1;  
                    idxRT = i;          %Store index of the RT such that the response time can be found
                    tau(numRT) = RTMTtau((find(RTMTtau(:,1) == idxRT)),3); %Find this index in RTMTtau matrix and get response time
                end
            end
        end
        ii = ii + 1;            %Increase an index by 1 for each new user
        getting_user_taus = [num2str(ii/24,'%.2f'),' %'] %There are 2400 individual users in present-data.
        tau = tau(tau ~= 0);    %Get rid of zeros from initialization.
        user_taus{1,ii} = u;    %ser(ii); %User-id
        user_taus{2,ii} = numTs;
        user_taus{3,ii} = numRT;
        user_taus{4,ii} = tau;  %Assign whole response time vector to cell
    end
% save('Results/TauDistr_IndividualUsers_detrended.mat'); %Important variables are user_taus
end