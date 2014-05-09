%================================
%=  Rustler Stefan, 2014        =
%=  <rustlers@student.ethz.ch>  =
%================================
%
% This simple function obtains the response times tau from a time series
% tevnts with the help of the matrix RTMTtau which contains the information
% on which time stamps are linked to give a tau. Since RTMTtau already
% contains taus in its third column. This function only makes sense to get
% taustars from tstars (since RTMTtau contains the ts as opposed to tstars)!
%
% Input:
% - RTMTtau_past: Matrix that contains the RT- and MT-indices w.r.t. the original tevnts-vector. It also contains the taus but we are interested in calculating the taus from a different set of tevnts (e.g. tstars).
% - tevnts: time series from which taus will be calculated. 
% 
% Output:
% - taus: response times that can, for example, be again plotted by using plotTaus(taus,'loglog','red')
% 
function [taus] = ts2taus(RTMTtau_past,tevnts)
    taus = zeros(100,1);
    for i = 1:(length(RTMTtau_past(:,1)))    %Loop through all MT-RT-links
        taus(i) = tevnts(RTMTtau_past(i,1)) - tevnts(RTMTtau_past(i,2));%Subtract the detrended timestamps of MT from RT to get a single tau
    end       
end