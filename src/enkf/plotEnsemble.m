function [varargout] = plotEnsemble(prm,xFull,c,fn,holdSP)
%% function [handle] = plotEnsemble(prm,xFull,c,fn,holdSP)
%
% Visualizes ensemble and ensemble mean of Q, S, G. 
%
% @param prm - system parametes (structure, see get_prmstruct.m)
% @param xFull - full states
% @param c - catchment identifier
% @param fn - figure number
% @param holdSP - Boolean hold subplots
% @return handle - handle to figure
%
% Usage:
%                 [handle] = plotEnsemble(prm,xFull,c,fn,holdSP)
%
% File:           plotEnsemble.m
%
% Created:        11/03/2013
%
% Last modified:  13/02/2013
%
% Author:         Tobias Siegfried (hydrosolutions ltd.)           
%
% Purpose:        Visualize ensemble mean and range for selected catchment.
%
% Description:    Main file of the data assimilation system.
%
% Revisions:      NA
%
% Copyright (C) 2013 Tobias Siegfried
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details.
%% startup
figure(fn)
catch2Show = c;

tStartShow = datenum(2013,04,25);

% PRECIPITATION AND POTENTIAL ET
pET = [prm.customprm.P(:,catch2Show) prm.customprm.ET(:,catch2Show)];
spH1 = subplot(5,1,1);
if holdSP 
    hold(spH1);
    plot(prm.forecasts.T,prm.forecasts.P,'b+-.');
    plot(prm.forecasts.T,prm.forecasts.ET,'g+');
    xlim([tStartShow prm.forecasts.T(end)]);
    datetick('x','dd/mm','keeplimits')
    %set(gca,'XTick',datenum(2013,01,01):1:prm.forecasts.T(end))
else
    plot(prm.customprm.time,pET) % here we cut off if GDAS is one day ahead of FEWS!
    legend('P','ET',2)
    xlim([min(prm.customprm.time) max(prm.customprm.time)])
    title(strcat('Precipitation and potential ET in Catchment', num2str(catch2Show)))
    ylabel('[mm]')
end

% RUNOFF
spH2 = subplot(5,1,2);
if holdSP 
    hold(spH2);
    xFullRS = reshape(xFull,[size(xFull,1) size(xFull,2) size(xFull,3) * size(xFull,4)]);
    minEnv = min(prm.customprm.qConvFac(c) * xFullRS(1 + (catch2Show - 1) * 4,:,:),[],3);
    maxEnv = max(prm.customprm.qConvFac(c) * xFullRS(1 + (catch2Show - 1) * 4,:,:),[],3);
    plot(prm.forecasts.T,mean(prm.customprm.qConvFac(c) * xFullRS(1 + (catch2Show - 1) * 4,:,:),3),'r-.')
    plot(prm.forecasts.T,minEnv,'k+-.')
    plot(prm.forecasts.T,maxEnv,'k+-.')
    xlim([tStartShow prm.forecasts.T(end)]);
    datetick('x','dd/mm','keeplimits')
    %set(spH2,'XTick',datenum(2013,01,01):2:prm.forecasts.T(end))
else
    hold on
    pol = [[prm.customprm.time;flipud(prm.customprm.time)] ...
    [min(prm.customprm.qConvFac(c) * xFull(1 + (catch2Show - 1) * 4,:,:),[],3)'; ...
    flipud(max(prm.customprm.qConvFac(c) * xFull(1 + (catch2Show - 1) * 4,:,:),[],3)')]];
    fill(pol(:,1),pol(:,2),[0.9 .9 .9])
    plot(prm.customprm.time,mean(prm.customprm.qConvFac(c) * xFull(1 + (catch2Show - 1) * 4,:,:),3),'r')
    hold off
    xlim([min(prm.customprm.time) max(prm.customprm.time)])
    legend('Ensemble spread','Q(mean)',2)
    title(strcat('Runoff in Catchment ', num2str(catch2Show)))
    ylabel('[m^{3}/s]')
end

% SOIL MOISTURE
spH3 = subplot(5,1,3);
if holdSP
    hold(spH3);
    minEnv = min(prm.customprm.qConvFac(c) * xFullRS(2 + (catch2Show - 1) * 4,:,:),[],3);
    maxEnv = max(prm.customprm.qConvFac(c) * xFullRS(2 + (catch2Show - 1) * 4,:,:),[],3);
    plot(prm.forecasts.T,mean(prm.customprm.qConvFac(c) * xFullRS(2 + (catch2Show - 1) * 4,:,:),3),'r-.')
    plot(prm.forecasts.T,minEnv,'k+-.')
    plot(prm.forecasts.T,maxEnv,'k+-.')
    xlim([tStartShow prm.forecasts.T(end)]);
    datetick('x','dd/mm','keeplimits')
    %set(spH3,'XTick',datenum(2013,01,01):2:prm.forecasts.T(end))
else
    hold on
    pol = [[prm.customprm.time;flipud(prm.customprm.time)] ...
    [min(xFull(2 + (catch2Show - 1) * 4,:,:),[],3)'; ...
    flipud(max(xFull(2 + (catch2Show - 1) * 4,:,:),[],3)')]];
    fill(pol(:,1),pol(:,2),[0.9 .9 .9])
    plot(prm.customprm.time,mean(xFull(2 + (catch2Show - 1) * 4,:,:),3),'r')
    hold off
    xlim([min(prm.customprm.time) max(prm.customprm.time)])
    legend('Ensemble spread','S(mean)',2)
    title(strcat('Soil Moisture in Catchment ', num2str(catch2Show)))
    ylabel('[mm]')
end

% GROUNDWATER
spH4 = subplot(5,1,4);
if holdSP 
    hold(spH4);
    minEnv = min(prm.customprm.qConvFac(c) * xFullRS(3 + (catch2Show - 1) * 4,:,:),[],3);
    maxEnv = max(prm.customprm.qConvFac(c) * xFullRS(3 + (catch2Show - 1) * 4,:,:),[],3);
    plot(prm.forecasts.T,mean(prm.customprm.qConvFac(c) * xFullRS(3 + (catch2Show - 1) * 4,:,:),3),'r-.')
    plot(prm.forecasts.T,minEnv,'k+-.')
    plot(prm.forecasts.T,maxEnv,'k+-.')
    xlim([tStartShow prm.forecasts.T(end)]);
    datetick('x','dd/mm','keeplimits')
    %set(spH4,'XTick',datenum(2013,01,01):2:prm.forecasts.T(end))
else
    hold on
    pol = [[prm.customprm.time;flipud(prm.customprm.time)] ...
    [min(xFull(3 + (catch2Show - 1) * 4,:,:),[],3)'; ...
    flipud(max(xFull(3 + (catch2Show - 1) * 4,:,:),[],3)')]];
    fill(pol(:,1),pol(:,2),[0.9 .9 .9])
    plot(prm.customprm.time,mean(xFull(3 + (catch2Show - 1) * 4,:,:),3),'r')
    hold off
    xlim([min(prm.customprm.time) max(prm.customprm.time)])
    legend('Ensemble spread','G(mean)',2)
    title(strcat('Groundwater in Catchment ', num2str(catch2Show)))
    ylabel('[mm]')
end

% ACTUAL ET
spH5 = subplot(5,1,5);
if holdSP 
    hold(spH5);
    minEnv = min(prm.customprm.qConvFac(c) * xFullRS(4 + (catch2Show - 1) * 4,:,:),[],3);
    maxEnv = max(prm.customprm.qConvFac(c) * xFullRS(4 + (catch2Show - 1) * 4,:,:),[],3);
    plot(prm.forecasts.T,mean(prm.customprm.qConvFac(c) * xFullRS(4 + (catch2Show - 1) * 4,:,:),3),'r-.')
    plot(prm.forecasts.T,minEnv,'k+-.')
    plot(prm.forecasts.T,maxEnv,'k+-.')
    xlim([tStartShow prm.forecasts.T(end)]);
    datetick('x','dd/mm','keeplimits')
    %set(spH5,'XTick',datenum(2013,01,01):2:prm.forecasts.T(end))
else
    hold on
    pol = [[prm.customprm.time;flipud(prm.customprm.time)] ...
    [min(xFull(4 + (catch2Show - 1) * 4,:,:),[],3)'; ...
    flipud(max(xFull(4 + (catch2Show - 1) * 4,:,:),[],3)')]];
    fill(pol(:,1),pol(:,2),[0.9 .9 .9])
    plot(prm.customprm.time,mean(xFull(4 + (catch2Show - 1) * 4,:,:),3),'r')
    hold off
    xlim([min(prm.customprm.time) max(prm.customprm.time)])
    legend('Ensemble spread','Etact(mean)',2)
    title(strcat('Actual ET in Catchment ', num2str(catch2Show)))
    ylabel('[mm/day]')
end

% RETURN HANDLES
varargout{1} = spH1; varargout{2} = spH2; varargout{3} = spH3; varargout{4} = spH4; varargout{5} = spH5;
end

