%================================
%=  Rustler Stefan, 2014        =
%=  <rustlers@student.ethz.ch>  =
%================================
% 
% Notes: 
% - This function plots the bottom subplots of the four graphs plotted by plot_RTs.m, i.e. RTMTtau(:,3) of that function is replaced by taus in this function.
% - taus always is RTMTtau(:,3), i.e. for both present and past data.
% - Unit of days was arbitrarily chosen. You can easily switch to hours by changing unit of taus and correcting axes labels.
% 
% Input:
% - taus: response times in DAYS, e.g. taus, taustars, AggregatedUserTausWithCutoff, or user_taus{4,uidx}, where uidx = 1...2400.
% - axes: either 'loglog' or 'semilog'
% - color: color of curve, e.g. 'blue' or 'green' or even [1 0.9 0.5] (relevant for functions that call this function several times for comparison, e.g. tstar2taustar.m
% 
% Output:
% - tau-histogram in units of days
% - tau-CCDF in units of days
%
function plot_taus(taus,axes,color)

    %Regular histogram on loglog-axis
    subplot(1,2,1);
        bins = 10000;
        g=max(taus)/bins*(1:bins);      %Rescale x-axis for histogram.
        h=hist(taus,bins);              %Need to specify binwidth! Also: We are performing 'absolute adding' (check algo-doc for more realistic options)
        loglog(g,h,'o','Color',color)   %Need to specify binwidth! Also: We are performing 'absolute adding' (check algo-doc for more realistic options)
        grid on;
        xlabel('ln \tau  (time between RT and MT in days)');
        ylabel('ln Count(X = \tau)');
        hold on;

    %Ranking plot on loglog- or semilogy-axis
    subplot(1,2,2);
        %y = taus;
        %y = sort(y,'descend');
        %x = (1:length(y)) / length(y); %Ranking plot (cdf, takes care of heavy-tail noise)
        [x,y] = ecdf(taus);
        x = 1-x;
        if strcmp(axes,'loglog')
            loglog(y,x,'.-','Color', color) %Note axes must be switched such that y-axis goes from 0 to 1.
            xlabel('ln \tau  (time between RT and MT in days)');
        end
        if strcmp(axes,'semilog')
            semilogy(y,x,'.-','Color',color) %Note axes must be switched such that y-axis goes from 0 to 1.
            xlabel('\tau  (time between RT and MT in days)');
        end
        ylabel('ln P(X \geq \tau)');
        %title([]);
        grid on;
end

