%================================
%=  Rustler Stefan, 2014        =
%=  <rustlers@student.ethz.ch>  =
%================================
%
% Input: 
% - user_taus:  Cell array returned by get_user_taus.m, containing user-ID, #tweets, #RTs and full vector of individual response times.
% - activity:   String of either 'tweets' or 'RTs' to state whether an "active" user is defined via number of RTs or tweets.
% - num:        Minimal number of tweets/RTs. NOTE: This can also be turned into an interval of #tweets (need to (un-)comment some parts in the code)
% - cutoff:     Number of last, i.e. longest, RTs (most extreme events) to be cut off when plotting and before aggregating the taus again. Note: the cutoff happens only for the users obtained with the criteria 'activity' and 'num'.
% - axes:       Can be either 'semilog' or 'loglog' for different plotting axes.
% 
% Output:
% - Plot with taus (response time) color-coded per user that fulfills the specified conditions
% - AggregUserTausWithCutoff: Same as user_taus only that "inactive" (as defined by 'activity' and 'num') users are disregarded as well as their extreme events, i.e. longest RTs, are cut off, as defined via 'cutoff'.
%
% Notes:
% - The following data contains user_taus and userIDs:
%    load('TauDistr_IndividualUsers.mat','user_taus')
%    load('TauDistr_IndividualUsers_detrended.mat','user_taus')
%

function [AggregUserTausWithCutoff] = plot_user_taus(user_taus,activity,num,cutoff,axes)
%% =====Find most active users=====

    %Find most active TWEETERS
    if strcmp(activity,'tweets')
        numTs_vec = [user_taus{2,:}];   %Vector that contains number of tweets for each user
        useridx = find(numTs_vec > num); %Vector that contains user index (NOT user-id! user-id = user_taus{1,useridx}) of users that tweeted more than 100 times.

    %     %Take the following four lines if you want to check LESS active users.
         useridx1 = find(numTs_vec > 0);
         useridx2 = find(numTs_vec < num);
         useridx = intersect(useridx1,useridx2); %Vector of idx that fulfill both conditions in the two previous lines
%          useridx = useridx(1:20);     %Since there are many of them useridx has to be cut down for the legend not to explode.    

    end

    %Find most active RETWEETERS
    if strcmp(activity,'RTs')
        numRT_vec = [user_taus{3,:}];
        useridx = find(numRT_vec > num);

    %     %Take the following four lines if you want to check LESS active users.
%          useridx1 = find(numRT_vec > 0);
%          useridx2 = find(numRT_vec < num);
%          useridx = intersect(useridx1,useridx2); %Vector of idx that fulfill both conditions in the two previous lines
%          useridx = useridx(1:20);     %Since there are many of them, useridx has to be cut down for the legend not to explode.  

    end

%% =====Plot results=====
%Uncomment line 54 to plot regular histogramm (comment some other lines,
%too, ofc.)

    AggregUserTausWithCutoff = [];
    i=0;
    legendcell = cell(1,50);    %Allocate some memory 
    h=figure;
    for uidx = useridx %Loop through users that fulfill condition given by num and activity
        i = i + 1;
        y = user_taus{4,uidx}; %This is basically the {taus} of user uidx
    %     h = hist(y,1000); loglog(h,'o','MarkerSize',13,'MarkerEdgeColor',[r g b]);
        y = sort(y,'descend');
        if cutoff < length(y)
            y = y((cutoff+1):length(y)); %Cutoff the last cutoff(#) events of each user (in useridx). Note: the last events correspond to the FIRST events in y!
        else
            y = [];
        end
        x = (1:length(y)) / length(y); %Ranking plot (cdf, takes care of heavy-tail noise)
    %     r = abs(sin(user_taus{1,uidx}/1000)); g = abs(cos(user_taus{1,uidx}/1000)); b = 0.5*abs(sin(user_taus{1,uidx}/1000)); %This only looks very complicated but it simply assigns to each userIDs a unique color!
        r = rand; g = rand; b = rand;

        if strcmp(axes,'semilog')
            semilogy(y,x,'.--','MarkerSize',13,'MarkerFaceColor',[r g b],'MarkerEdgeColor',[r g b]) %Note axes must be switched such that y-axis goes from 0 to 1.
            xlabel('\tau  (time between RT and MT in hours)');
        end
        if strcmp(axes,'loglog')
            loglog(y,x,'.--','MarkerSize',13,'MarkerFaceColor',[r g b],'MarkerEdgeColor',[r g b]) %Note axes must be switched such that y-axis goes from 0 to 1.
            xlabel('ln \tau  (time between RT and MT in hours)');
        end
        
        legendcell{1,i} = ['User-ID: ',num2str(user_taus{1,uidx}), ' (#Tweets = ',num2str(user_taus{2,uidx}),', #RTs =',num2str(user_taus{3,uidx}),')'];
        hold on;
        
        AggregUserTausWithCutoff = horzcat(AggregUserTausWithCutoff,y); 
        
    end 
    
    legendcell(cellfun(@isempty,legendcell)) = []; %Delete empty cell entries because they are of type non-string and not handable by legend()!
    legend(legendcell);
    ylabel('ln P(X \geq \tau)');
    grid on;
    hold off;
    
    saveas(h,[pwd,'/Results/Detrended_user_statistics_of_present.pdf']);
end