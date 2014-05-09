function [timeperunit, count] = get_timeperunit(tevnts,unit)

    %INPUT:
    %- tevnts: time stamps to be split into 'unit's (both pres and past data!). Format needed: 'yyyy-mm-dd HH:MM:SS'. Needs to be prepared by previos function!
    %- unit: Can be 'month' or 'day'. week is NOT implemented yet!
    %
    %OUTPUT:
    %- timeperunit: cell array with each entry being the {tevnts} of the unit, i.e. day, week or month
    %- count: number of events in the whole series to check whether they are conserved.
    %
    %EXAMPLE:
    %[timeperunit] = get_timeperunit_present(tALL,'month');

    %% ===== Load and write data =====

    SetSize = size(tevnts,1);
    date = zeros(SetSize,1); %Preallocate memory

    for i=1:SetSize
        date(i) = datenum(tevnts(i),'yyyy-mm-dd HH:MM:SS');
%         t = datevec(presC{1}(i),'yyyy-mm-dd HH:MM:SS'); %This is just an auxilliary
%         variable (only needed in the next step)
%         time(i) = datenum(horzcat([1988 6 9],t(4:6))); %Set same date for all since we only want clock time
    end

    %Put 1dim vector date into vectors for each unit (timeperunit)
    timeperunit = cell(30,1); %Preallocate memory
    j = 1;
    count = 0;
    
    if strcmp(unit,'day')
        day = 1;
        for i=1:(length(date)-1)
            a = datevec(date(i+1));
            b = datevec(date(i));
            if  a(3) > b(3) || a(2) > b(2) % The second OR-condition accounts for the month border!
                timeperunit{day} = date(j:i);
                j = i + 1;
                b(4:6) = [23 59 59]; %Set i-th event to just before i+1-th event s.t. ceil() is always correct. (Without this there are cases where one would need floor() instead, e.g. when last tweet of the unit is very early.)
                day = day + ceil(date(i+1)-datenum(b)); %Add the unit difference (thus floor to get round units, NOT exact difference) because it can happen that units are skipped, i.e. there are no tweets!
%                 tag(i) = day;
%                 if ceil(date(i+1)-date(i)) > 1 %This is a different count: number
%                 of days skipped.
%                     count = count + 1;
%                 end
            end
        end
        timeperunit{day} = date(j:length(date)); % Manually fill up last incomplete unit as it is not captured by the code snippet just before.
        
%         %NOTE: #units=936 #times a unit was skipped = 76 (count) in this code even though
%         %datenum(presC{3}(SetSize))-datenum(presC{3}(1))=939.16 (->941). So 936+76 > 941, which is weird. All in all,
%         %though, everyunit is captured it might only be that their enumeration
%         %is false, i.e. #skippedunits is too large -> too many empty-valued
%         %units.
        
        for i=1:day
            count = count + numel(timeperunit{i});
        end
    end
    
    if strcmp(unit,'month')
        month = 1;
        for i=1:(length(date)-1)
            a = datevec(date(i+1));
            b = datevec(date(i));
            if  a(2) > b(2) || a(1) > b(1) % The second OR-condition accounts for the month border!
                timeperunit{month} = date(j:i);
                j = i + 1;
                month = month + 1;
            end
        end
        timeperunit{month} = date(j:length(date));
        for i=1:month
            count = count + numel(timeperunit{i});
        end
    end
%     NumberOfEvents = length(date);    %Count should agree with this value (#events). It does for 'month' but not for 'day' for some reason.
end