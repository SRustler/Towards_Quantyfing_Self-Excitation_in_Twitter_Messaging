%================================
%=  Rustler Stefan, 2013        =
%=  <rustlers@student.ethz.ch>  =
%================================
%
% NOTEs: 
% - To speed up tests with this function (only relevant when visuals=1), the lines for the Hawkes process under subplot(4) can be commented out. 
% - Smoothing parameters have to be changed within the function, i.e. they are not taken from the function arguments.
% - When plotting data, the figures will be stored in to /Results... . Thus Make sure you are in the parent directory /Code.
% - The lines commented out and including pTmax are obselete after having introduced firstmom and firstevnt. They are not removed yet in case we need to go back.
%
% Input:
% - binsize: always in minutes even if you choose to loop through months or weeks instead of days. E.g. 1,10,60.
% - unit: currently only 'day','month'. weekly not yet implemented.
% - idx: index for ii-loop, i.e. through which units, e.g. days, you want to loop. 
% - tMT: time stamps of MT-events (RT-filtered). Format needed: 'yyyy-mm-dd HH:MM:SS'.
% - tALL: time stamps of all events (unaltered). Format needed: 'yyyy-mm-dd HH:MM:SS'.
% - visual: 1 if you want each unit to be visualized (no figures will be saved), 0 if you are just interested in tstartot and gofcell
% 
% Output:
% - 4Graphs per unit: R(t) with spline fit, R(t*) with mean of spline fit, N(t), residuals with quantiles
% - tstartot: total transformed time t* in units of hours (Note: date returned by ActVsTime_present.m is in units of days!)
% - gofcell: cell array of struct of gof parameters of the spline fit.
%
% ISSUEs: 
% - Minor: Sometimes tstar entries at the end are NaN. This sometimes creates cluster of NaN's, i.e. NaN-gaps in the N(t)-curve. WHY??
% - Verify code with visuals=1 after having introduced firstmom and firstevnt.

function [tstartot,gofcell] = detrend_anyFit(binsize,unit,idx,tMT,tALL,visual)%,tMTs_perunit,tALL_perunit)
    %% ===== Load data and declare basic variables =====
    
    [tALL_perunit,~] = get_timeperunit(tALL,unit);  %Now tevtns per unit can be called by timeperunit{i} (currently only days or months)(Unit: #days since reference 1-Jan-00)
    [tMTs_perunit,~] = get_timeperunit(tMT,unit);   %Now tevtns per unit can be called by timeperunit{i} (currently only days or months)(Unit: #days since reference 1-Jan-00)    
        
    if strcmp(unit,'day')
        Tmax = 24;              %Hours of the day
%         pTmax = Tmax;           %Tmax of previous day/month. (We need to let the t* of the next day/month start at pTmax!)
        smoothingParam = 0.3;   %Smoothness parameter for spline fitting. 0 ist ganz grob, 1 ganz fein.
        t = 0:(1/60*binsize):(Tmax-(1/60*binsize)); %Change units from minutes = [binsize] to hours!
    end

    month = false;  %default state
    if strcmp(unit,'month')
        month = true;
        smoothingParam = 0.8;   %Smoothness parameter for spline fitting. 0 ist ganz grob, 1 ganz fein. 
    end
    
    firstmom = datevec(tALL_perunit{1}(1));      
    firstmom(4:6) = [0 0 0]; firstmom = 24*datenum(firstmom);  %
%     lastmom = firstmom;
    tstartot = zeros(10000,1); %Initialize vector of all t*.
    gofcell = cell(numel(idx),1); %Initialize cell with gofs for each detrended unit.
    
    ['Time stamps split into ',unit,'s. Now detrending ...'] %Print this to check progress of code.

    %% ===== Detrending ====
    if visual == 0
        for ii = idx
            
            if numel(tALL_perunit{ii}) == 0 || numel(tMTs_perunit{ii}) == 0
                continue    %In present data it can happen that some days do not have any tweets. Continue to next day/month in this case.   
            end
            
%             lastmom = datevec(tALL_perunit{ii}(1));        %Get datevec of previous unit.
%             lastmom(4:6) = [0 0 0]; lastmom = datenum(lastmom); %Get last moment of previous month.
            
            if month   %If unit='month', Tmax and pTmax vary for each ii!
                Y = datevec(tALL_perunit{ii}(1)); Y = Y(1); %Get current year for time-axis normalization: Check 1st time stamp of current month and get the year of it, i.e. 1st entry of datestring.
                M = datevec(tALL_perunit{ii}(1)); M = M(2); %Get current month for time-axis normalization: Check 1st time stamp of current month and get the year of it, i.e. 2nd entry of datestring.
%                 pTmax = 0;
%                  if ii>1 
%                      pY = datevec(tALL_perunit{ii-1}(1)); pY = pY(1); %Y of previous month: needed in line before concatenation of tstar
%                      pM = datevec(tALL_perunit{ii-1}(2)); pM = pM(2); %M of previous month: needed in line before concatenation of tstar
%                      pTmax = eomday(pY,pM) * 24;                      %Get last day of month pM/pY. Multiply by 24 to get number of hours in that month
%                      lastmom = datevec(tALL_perunit{ii-1}(1));        %Get datevec of previous unit.
%                      lastmom(4:6) = [23 59 59]; lastmom = datenum(lastmom); %Get last moment of previous month.
%                  end
                Tmax = eomday(Y,M) * 24;
                t = 0:(1/60*binsize):(Tmax-(1/60*binsize)); %We need this line again for the case of months.
            end
            
            %Detrending:
            tevnts = tMTs_perunit{ii};          %Load only MTs for fitting purposes of {t_i}_unit. Recall: mu(t) as unconditional intensity does not depend on endo-events, i.e. RTs!
            tevnts = (tevnts - tevnts(1)) * 24; %Now match scale to [0;Tmax[ HOURS, since datenum's are ALWAYS in days.
            h = hist(tevnts,numel(t));          %Change {t_i}_unit to R(t)_unit with correct bin size. Thus hist(data,nbins).
            [splinefn,gof,~] = fit(t.',h.','smoothingspline','smoothingParam',smoothingParam); %Perform spline fit.
            int = integrate(splinefn,t.',0);    %Basically assign to each bin the integration-value until this bin. So numel(int)=numel(t).
            K = max(int) / max(tevnts);         %Normalization factor
            tevnts = tALL_perunit{ii};          %Load ALL tweets for detrending!
            firstevnt = (tevnts(1)) * 24;       %Store time of first event in unit 
            tevnts = 24*tevnts - firstevnt;     %Now match scale to [0;24[ (or e.g. [0;24*31[ for months) since datenum's are ALWAYS in days.
            xi = tevnts;                        %Detrend tevnts of ALL tweets!
            yi = interp1(t,int,xi);%,'cubic');  %int has viewer points than numel(t). Interpolate such that point of tevnts that would fall between points of t get their correct respective y-value!
            tstar = yi/K;                       %Normalize t*. Now we got one tstar for each tevnts of tALL!
            gofcell{ii} = gof;                  %Store gof-value for this particular day/month.
%             tstar = tstar + (ii-1)*pTmax;       %Since it was set that 1=1h and then normalized to one day (i.e. Tmax of the PREVIOUS month = pTmax = eomday(Y,M-1)*24), 24h (Tmax s.t. also 'month' works) have to be added manually such that tstartot is monotonically increasing.            
            tstar = tstar + firstevnt-firstmom; %Calibrate w.r.t first moment of current unit:
            
            tstartot = vertcat(tstartot,tstar); %QUESTION: Why is the order of magnitude for tstar higher (780)?
        end
    end

%% ===== Detrending with plotting ====

    if visual == 1
        for ii=idx
            
%=========CALENDAR VIEW====start
%Uncomment the following and comment the rest in for ii=idx to plot in
%calendar view
%             subplot(6,7,ii+4); %Plot unconditional intensity in calendar view. Add 4 in order to start at 5 because day 1 happens to be a Friday. This lets the week start on a Monday.
%                 h = hist(tALL_perunit{ii},numel(t)); %h = hist(sort(ttime),numel(t));
%                 bar(t,h);hold on;   %Plot ALL tweets
%                 tevnts = tMTs_perunit{ii};    %Load only MTs for fitting purposes. Recall: mu depends on exo-events, i.e. MTs, only!
%                     [~,DayName] = weekday(tevnts(1));   %Extract weekday
%                     a = datestr(tevnts(1),'dd-mm-yy');  %Extract date
%                     tevnts = (tevnts - tevnts(1)) * 24; %Now match scale to [0;24[
%                 h = hist(tevnts,numel(t)); %h = hist(sort(ttime),numel(t)); %Perform fit only on MTs!
%                 [splinefn,gof,~] = fit(t.',h.','smoothingspline','smoothingParam',smoothingParam); %Perform fir 
%                 plot(splinefn,'r-');
%                 ylim([0 max(h)+10]) %Add five such that it does not look it's cut off
%                 xlim([0 Tmax]); grid on;
% 
%                 title(['Day ', num2str(ii),' (',DayName,', ',a,')']);
% %             
%=========CALENDAR VIEW====end
           
            if numel(tALL_perunit{ii}) == 0 || numel(tMTs_perunit{ii}) == 0
                continue    %If unit='day' it can happen that some days do not have any tweets.   
            end
            
            
            if month        %If unit='month', Tmax varies for each ii!
                Y = datevec(tALL_perunit{ii}(1)); Y = Y(1);
                M = datevec(tALL_perunit{ii}(2)); M = M(2);
%                 pTmax = 0;
%                 if ii>1
%                     pY = datevec(tALL_perunit{ii-1}(1)); pY = pY(1); %Y of previous month: needed in line before concatenation of tstar
%                     pM = datevec(tALL_perunit{ii-1}(2)); pM = pM(2); %M of previous month: needed in line before concatenation of tstar
%                     pTmax = eomday(pY,pM) * 24;
%                 end
                Tmax = eomday(Y,M) * 24;
                t = 0:(1/60*binsize):(Tmax-(1/60*binsize)); %Change units from minutes = [binsize] to hours!
            end
         
            
            scrsz = get(0,'ScreenSize');
            fig1 = figure('Position',[1 scrsz(4)/1.4 scrsz(3)/1.4 scrsz(4)/1.4]); %Adapts figure to screen size.
            set(fig1,'PaperPositionMode','auto'); %Saves figure as displayed
            
            
            subplot(2,2,1); %Plot fitted unconditional intensity on R(t)
                H = hist(tALL_perunit{ii},numel(t));    %h = hist(sort(ttime),numel(t));
                h = hist(tMTs_perunit{ii},numel(t));    
                bar(t,[H' (H-h)'],'stack');hold on;     %Plot MTs and separately stack the RTs on top of each bar.
                tevnts = tMTs_perunit{ii};              %Load only MTs for fitting purposes. Recall: mu depends on exo-events, i.e. MTs, only!
                    [~,DayName] = weekday(tevnts(1));   %Extract weekday
                    a = datestr(tevnts(1),'dd-mm-yy');  %Extract date
                    tevnts = (tevnts - tevnts(1)) * 24; %Now match scale to [0;24[
                h = hist(tevnts,numel(t)); %h = hist(sort(ttime),numel(t)); %Perform fit only on MTs!
                [splinefn,gof,~] = fit(t.',h.','smoothingspline','smoothingParam',smoothingParam); %Perform fir 
                plot(splinefn,'r-');
                hold off;
                
                ylim([0 max(h)+10]) %Add five such that it does not look it's cut off
                xlim([0 Tmax]); grid on;
                title(['Day ', num2str(ii),' (',DayName,', ',a,')']);
                xlabel(['Day time t [h in',num2str(binsize),'-min-bins]']);
                ylabel(['Event rate R(t) [tweets per ',num2str(binsize),' min]']);
                legend('RT-filtered process R_{MT}(t)','Stacked RTs','Spline fit to \mu(t)','Location','NorthWest');   
                if strcmp(unit,'month') %If true, previous labels will simply be over-written.
                    title(['Month ', num2str(ii),' (First day: ',a,')']);
                    xlabel(['Time t [h in ',num2str(binsize),'-min-bins]']);
                end


            subplot(2,2,3); %Time transform
                int = integrate(splinefn,t.',0);
                K = max(int) / max(tevnts);     %Normalization factor
                tevnts = tALL_perunit{ii};      %Load ALL tweets for detrending!
                firstevnt = (tevnts(1)) * 24;   %Store time of first event in unit 
                tevnts = 24*tevnts - firstevnt; %Now match scale to [0;24[ (or e.g. [0;24*31[ for months) since datenum's are ALWAYS in days.                xi = tevnts;                %Detrend tevnts of ALL tweets!
                xi = tevnts;                    %Detrend tevnts of ALL tweets!                
                yi = interp1(t,int,xi);%,'cubic'); %Interpolate such that point of tevnts that would fall between points of t get their correct respective y-value!
                tstar = yi/K;                   %Normalize t*
                hstar = histc(tstar,t);
                bar(t,hstar);hold on            %Transformed rate
                mu = feval(splinefn,t.');       %Convert fitobject to vector whose mean can be calculated
                plot(t,mean(mu),'g.-');         %Calculate and plot mean
                hold off;
                
                ylim([0 max(h)+10])         %Same limit as for original time for better comparability
                xlim([0 Tmax]); grid on;
                xlabel(['Transformed day time t* [',num2str(binsize),' min]']);
                ylabel(['Transformed event rate R(t*) [tweets per ',num2str(binsize),' min]']);
                legend('Process with transformed time \lambda(t*)',['Unconditional constant intensity \mu = ',num2str(mean(mu))],'Location','NorthWest');
                if strcmp(unit,'month') %If true, previous labels will simply be over-written.
                    xlabel(['Transformed time t* [h of',num2str(binsize),'-min-bins]']);
                end                
                
                
            subplot(2,2,2); %Cumulative number of events
                y = 1:length(tstar);
                n = max(tstar)/max(y);  %Normalization such that y-axis matches x-axis in scale. Recall that x-axis (tstar) was normalized before, too!
                y = y * n;
                plot(tevnts,y,'b-')     %Original complete process. Note that this corresponds exactly to F(t) which is the primitive function of lambda(t). Test: %plot(t, F(t),'r', 'linewidth', 2)
                hold on;
                plot(t,24/max(int)*int,'r-'); %Original "unconditional" (spline-fit) process. This corresponds to the primitive function of the assumed mu(t).
                plot(tstar,y,'g-')      %Transformed process
                plot(y,y,'--k')         %Perfect Poisson process with lambda=1
                hold off; 
                
                legend('Original complete process with \lambda(t)','Original ``unconditional'' process with \mu(t)','Transformed complete process with \lambda(t*)',['Poisson process with \lambda =',num2str(1/n)],'Location','NorthWest') %NOTE: We dont get a uniform process with lambda=1 because of enforced normalization with K!
                xlim([0 Tmax]); grid on;
                xlabel('Time [h]')
                ylabel(['Counting process N [1/',num2str(1/n),']']) 

                
            subplot(2,2,4);
                %The following operation takes very long!
                H=hawkes('exp');
                H.plot_residuals(tstar(tstar<Tmax));    %Sometimes it happens that tstar has some NaN entries. This will result in an error at this point.
                xlim([0 Tmax]); grid on;

                
            gofcell{ii} = gof;
            tstar = tstar + firstevnt-firstmom; %Calibrate w.r.t first moment of current unit:
            
%             %ISSUE: This if-part is a cheat! It corrects some unidentified flaw of the previous line, which somehow does not always make tstartot monotonically increasing.
%             if tstartot(end)>tstar(1)   
%                 tstar = tstar + (tstartot(end)-tstar(1));
%             end
%             %Flaw seems to be resolved! Only issue left: sometimes there NaN-gaps between units.
            
            tstartot = vertcat(tstartot,tstar);  %QUESTION: Why is the order of magnitude for tstar higher (780)?
            saveas(fig1,['Results/Rates_',unit,num2str(ii),'.eps'],'epsc')
        end
    end
    
    tstartot=tstartot((find(tstartot>0)):end); %Cut off first non-empty value from initialization as well as zeros in the beginning (Why did they suddenly come?)
    whos tstartot %Check whether Size agrees with number of tweets in {t_i} (Recall: Detrending is a 1:1-mapping!).
end