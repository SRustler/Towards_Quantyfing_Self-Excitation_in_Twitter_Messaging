%================================
%=  Rustler Stefan, 2014        =
%=  <rustlers@student.ethz.ch>  =
%================================
%
% Issues:
% - Retweets that do not follow the standard way of retweeting, "... RT @RTd_user: MotherTweet ..." (or space instead :), will not be correctly handled by this function!
% - In particular if 'RT @user' is put at the end of the RT, the reconstructed MT will be wrong.
% - Space completion for usernames and tweets happens for each Levenshtein distance, which is not wrong but inefficient.
%
% Data snippet:
%   pastC{1}(i) - time stamp|pastC{2}(i) - tweet               |pastC{3}(i) - username
%   '2009-10-23 13:22:07'   |'RT @koprax: Bin stark...         |'michelreimon'        
%   '2009-10-23 13:23:37'   |'Werde mich gleich ... #unibrennt'|'staenkerliese'       
%   '2009-10-23 13:24:47'   |'Wir demostrieren ... #unsereuni' |'koprax'              
%
% Input:
% - pastC: Cell of past data that contains time stamp, tweet message, username (%s %s %s).
% - LevenshteinDist: Vector of Levenshtein distances to loop through. Minimum distance is zero (no agreement of strings), the maximum 140 (full agreement of 140-character tweets). The distance defines the criterion for a potential MT to be accepted as actual MT of a RT found.
%
% Output:
% - RTMTtau_pastvsLeven: Cell of past data that for each LevenshteinDist-element contains the matrix with RT-index, the respective MT-index and response times between the MT and its RT (in units of days). 
%
% Notes:
% - In case there are multiple 'RT @' in the tweet, the right most one will be taken, assuming that this refers to the most original MTer. However, there are some case where the 'RT @RTd_user' is put at the end!
% - RTd_user will be the last successfully determined re-tweeted user.
% - user(MT) will be the the re-tweeter that was just found.
% - Cell syntax:
%   pastC{3}(1) = 'Peter' => cell with numel=1
%   pastC{3}{1} = Peter   => string with numel=5

function [RTMTtau_pastvsLeven] = get_taus_past(pastC,LevenshteinDist)

    RTMTtau_past=zeros(floor(0.25*size(pastC{1},1)),3); %Reserve some space for matrix with RT-index, MT-index and response times between the MT and its RT.
    RTMTtau_pastvsLeven = cell(1,numel(LevenshteinDist)); %Reserve space for cell containing RTMTtau_past for each LevenshteinDist.
    
    levidx = 0;                         %Levenshtein index to replace calculations. E.g. leven = [0:20:140] --> levidx = [1:1:8].
    for leven = LevenshteinDist 
        ii = 0;                         %RT-idx starting at zero and ending at total number of RTs found.
        for i=1:size(pastC{1},1)        %Loop through all tweets to search for the character sequence 'RT @'.
            diff = 15 - numel(pastC{3}{i});             %Ensure that all the USERNAMES are of same string length.
            if diff ~= 0
                diff = char(zeros(1,diff));             %Create some string of length diff.
                pastC{3}{i} = [pastC{3}{i} diff];       %Complement USERNAME to fill maximal number of chars with spaces.
            end
            
            diff = 141 - numel(pastC{2}{i});        %Ensure that all the TWEETS are of same string length.
            if diff ~= 0
                diff = char(zeros(1,diff));         %Create string of length diff.
                pastC{2}{i} = [pastC{2}{i} diff];   %Complement TWEET to fill maximal number of chars with spaces.
            end
            %ISSUE: We are doing this completion for all leven distances. This is inefficient!
            
            tweet_str = pastC{2}{i};  
            b = max(strfind(tweet_str,'RT @'));     %Get index of position of where 'RT @' occurs in the tweet last (thus 'max'). This accounts for the possibility to have several 'RT @' in the tweet by choosing the most originial mother tweeter.

            %% Get retweeted user (RTd-user)
            if numel(b)~=0 %Get username of RTd-user (=mother tweeter) if there is at least one 'RT @' in the tweet.
                
                %Strategy to find the RTd user exploits the fact that RTd-user occurs after 'RT @' and before the first ' ' or ':' after that.
                time_RT = datenum(pastC{1}(i),'yyyy-mm-dd HH:MM:SS');   %Record time of RT for response time calculation. datenum(string, 'strformat') saves date in numerical value in units of days (since 1.1.1900).
                RTd_user = [tweet_str((b+4):numel(tweet_str)),' '];     %First get rid of 'RT @' and everything before that. Then artificially add ' ' to account for 'RT @RTd-user' occurring at the end of the tweet. RTd-user at this point is to be splitted into the actual RTd-user and the MT.

                %Find end-index of RTd-user:
                b1 = strfind(RTd_user,' '); %Get indices of position of ' '
                b2 = strfind(RTd_user,':'); %Get indices of position of ':'
                if numel(b2)==0             %The tweet does not always contain ':'!
                    b2 = b1(1) + 1;         %Make sure b2 is bigger and therefore not chosen for b3.
                end    
                b3 = max(b1(1),b2(1));      %Get index+1 of position of where the RTd username ends. 'max' makes sure that the righ most, i.e. most probably the most original, MTer is chosen.

                MT = RTd_user((b3+1):(numel(RTd_user)-1));  %Obtain MT and subtract 1 to discard ' ' that was added before.
                %ISSUE: If 'RT @RTd-user' really occurs at the end of the tweet sometimes, then MT cannot be right in this case.

                
                diff = 140 - numel(MT);         %Ensure that all the RECONSTRUCTED MTs are also of same string length.
                if diff > 0
                    diff = char(zeros(1,diff)); %Create string of length diff
                    MT = [MT diff];             %Complement MT to fill maximal number of chars with spaces
                end
                if diff < 0
                    MT = MT(1:141); %In case MT is longer than 140 (somehow that happened!). 
                end 

                RTd_user = RTd_user(1:(b3-1));  %Secondly get rid of everything after the RTd username.

                diff = 15 - numel(RTd_user);    %Make sure the RTd_user is also of length 15
                if diff > 0
                    diff = char(zeros(1,diff)); %Create string of length diff
                    RTd_user = [RTd_user diff]; %Complement USERNAME to fill maximal number of chars with spaces
                end
                %ISSUE: RTd_user at this point is by mistake sometimes longer than 15 chars in which case the following line applies:
                if diff < 0
                    RTd_user = RTd_user(1:15);  %In case RTd_user is longer than 15.
                end   
                %RTd_user = [RTd_user char(zeros(1,(15 - numel(RTd_user))))]; %Finally fill up with char 'space' s.t. that maximum length is
                %fulfilled
                
                %% Find original MT and get time of MT
                for j=i-1:-1:1 %Scan whole history of tweets previous to RT, starting from most recent in history
                    if floor(sum(pastC{3}{j}==RTd_user)/15) && strdist(pastC{2}{j},MT) <= leven
%                     if strcmp(pastC{3}{j},RTd_user) && strdist(pastC{2}{j},MT) <= leven %&& numel(MT)==numel(pastC{2}{j}) && floor(sum(pastC{2}{j}==MT)/numel(MT))%First make sure the numel is equal (otherwise you cannot compare strings). Then compare string entry by entry. If strings are equal, 1 occurs for each equal char. Sum up, divide by number of entries to get exactly 1. If only one char is NOT equal, the floor function will return ZERO.

                       time_MT = datenum(pastC{1}(j),'yyyy-mm-dd HH:MM:SS'); %I.e.: Check if MT is contained in pastC{3}{k}

                       %Do the following only if a MT was found for a RT
                       ii = ii + 1;                             %RT-counter
                       RTMTtau_past(ii,1) = i;                  %Assign index of RT w.r.t. originial set (pastC)
                       RTMTtau_past(ii,2) = j;                  %Assign index of MT w.r.t. originial set (pastC)
                       RTMTtau_past(ii,3) = time_RT - time_MT;  %Response time between MT and RT. This in unit of days!

%                       [RTd_user,'`s tweet: "',MT,'" was retweeted at ',num2str(time_RT)]; %Check RTd_user and their MT.

                       break; %Break k-for-loop once tau has been obtained and thus prevent time_MT to be overwritten by earlier message of RTd_user.
                    end     
                end
            end 
            getting_RTs_for_LevenshteinDist = [num2str(i/511.14,'%.3f'),'%    | leven = ',num2str(leven)] %Print progress of code.
        end 
        RTMTtau_past = RTMTtau_past(1:ii,1:3); %Discard unoverwritten zeros such that they are not plotted in histogram!
        levidx = levidx + 1;
        RTMTtau_pastvsLeven{1,levidx} = RTMTtau_past; 
%         save([pwd,'/Results/RTs_vs_Levenshtein/Data_LevenshteinDist',num2str(leven),'.mat'],'RTMTtau_past')
    end
end