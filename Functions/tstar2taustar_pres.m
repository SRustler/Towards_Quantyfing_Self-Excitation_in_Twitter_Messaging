% This function compares the detrended and original tau distribution.
%
% Input:
% - tstartot: detrended time series of events
% - str: either 'loglog' or 'semilog'
% 
% Output:
% - taustars: transformed response times that can, for example, be again plotted by using plotTaus(taus,'loglog','red')
% 
% Example:
% [taus] = tstar2taustar([30 120], tstartot, 'loglog');

function [taustars] = tstar2taustar_pres(tstartot,RTvsDt,str)

    taustars = zeros(1000,1);
    a = RTvsDt(:,1);    %Vector of RT indices (per row larger than the MT index)
    a=a(a<=numel(tstartot));    %tstartot might have not as many elements as there are events (e.g. in present-data). Do not consider the taus for indices that exceed numel(tstartot).
    for i = 1:(length(a)) %Loop through all MT-RT-links
        taustars(i) = 1/24 * (tstartot(RTvsDt(i,1)) - tstartot(RTvsDt(i,2))); %Subtract the detrended timestamps of MT from RT to get a single tau. Divide by 24 to get units of days.
    end
%     figure
    plot_taus(1/24*RTvsDt(:,3),str,'red'); %Original taus in units of hours(!)
    hold on;
    plot_taus(taustars,str,'green'); %Transformed taus in units of days(!)
    legend('Original time \tau','Transformed time \tau*','Location','SouthWest')
    clear RTvsDt
   
    hold off
end