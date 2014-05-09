%================================
%=  Rustler Stefan, 2014        =
%=  <rustlers@student.ethz.ch>  =
%================================
%
% Input:
% - RTMTtau is the cell that contains per row: idxRT, idxMT, tau(responsetime)
% - SetSize is the size of the original data, i.e. # all tweets (51,114 for past). Note that size(RTMTtau) < SetSize!
% - str: choose between 'loglog' or 'semilog'
%
% Output:
% - Prominent four plots visualizing RTMTtau: idxMT vs tau, idxMT vs idxRT, tau-histogram, tau-CCDF. 
%
function plot_RTs(RTMTtau,SetSize,str)

    subplot(2,2,1);plot(RTMTtau(:,3),RTMTtau(:,2),'o')
        xlim([0 Inf]);
        xlabel('Time between RT and MT in days');
        ylabel('Index of MT'); ylim([0 SetSize]); 
        grid on;

    subplot(2,2,2);plot(RTMTtau(:,1),RTMTtau(:,2),'o')
        xlabel('Index of RT'); xlim([0 SetSize]);
        ylabel('Index of MT'); ylim([0 SetSize]); 
        grid on;

        %subplot(2,2,3);hist(RTMTtau(:,3),100);
        %xlabel('Time between RT and MT in days');
        %ylabel('Count of times until RT'); 
        %grid on;

        
    subplot(2,2,3); 
        bins = 10000;
        g=max(RTMTtau(:,3))/bins*(1:bins);    %Rescale x-axis for histogram.
        h=hist(RTMTtau(:,3),bins); %Need to specify binwidth! Also: we are performing 'absolute adding' (check algo-doc for more realistic options)
    
        loglog(g,h,'mo') %Need to specify binwidth! Also: we are performing 'absolute adding' (check algo-doc for more realistic options)
        grid on;
        xlabel('ln \tau  (time between RT and MT in days)');
        ylabel('ln Count(X = \tau)');

    subplot(2,2,4);
        %bins = 10000;
        %g=max(RTMTtau(:,3))/bins*(1:bins);    %Rescale x-axis for histogram.
        %h=hist(RTMTtau(:,3),bins); %Need to specify binwidth! Also: we are performing 'absolute adding' (check algo-doc for more realistic options)
        %loglog(g,h,'x')%hist(RTMTtau(:,3),100); %Need to specify binwidth! Also: we are performing 'absolute adding' (check algo-doc for more realistic options)
        %loglog(g,h,'x')%hist(RTMTtau(:,3),100); %Need to specify binwidth! Also: we are performing 'absolute adding' (check algo-doc for more realistic options)
        y = RTMTtau(:,3);
        y = sort(y,'descend');
        x = (1:length(y)) / length(y); %Ranking plot (cdf, takes care of heavy-tail noise)
        if strcmp(str,'loglog')
            loglog(y,x,'m.-') %Note axes must be switched such that y-axis goes from 0 to 1.
            xlabel('ln \tau  (time between RT and MT in days)');
        end
        if strcmp(str,'semilog')
            semilogy(y,x,'m.-') %Note axes must be switched such that y-axis goes from 0 to 1.
            xlabel('\tau  (time between RT and MT in days)');
        end
        ylabel('ln P(X \geq \tau)');
        grid on;

end

