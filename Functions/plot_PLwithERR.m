function [L] = plot_PLwithERR(RTvsDtCOL3)
%x = RTvsDt(:,3); %In plotRTs.m this is y
x = RTvsDtCOL3;

[alpha, xmin, L] = plfit(x);    %Get get best alpha, lowest xmin, and corresponding likelihood
[dalpha, dxmin, n] = plvar(x,'xmin',xmin); %Get uncertainties for previously determined xmin
dalpha = roundsd(dalpha,2); %Two significant digits for error
dxmin = roundsd(dxmin,2);  %Two significant digits for error

plplot(x, xmin, alpha, dxmin, dalpha);

end
