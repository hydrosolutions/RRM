%% script ASSIMILATION TRAINING
%
% This script performs assimilation with training data from NCEP CFS 
% for the Themi.
%
% @param fName      - Loads basin-specific training data
% @param dataPath   - Path where basin-specific training data is stored.
% @param savePath   - Output path
% @param saveName   - Output name
%
% Usage:
%                 NA
%
% File:           startAssimilationTraining.m
%
% Created:        01/11/2012
%
% Last modified:  08/05/2013
%                   - Training restored.
%
% Author:         Tobias Siegfried (hydrosolutions ltd.)           
%
% Purpose:        Model training script
%
% Description:    
%
% Trains a 'blank' RR model with NCEP CFS data (see: Saha, Suranjana, and
% Coauthors, 2010: The NCEP Climate Forecast System Reanalysis. Bull. Amer.
% Meteor. Soc., 91, 1015.1057. doi: 10.1175/2010BAMS3001.1)
%
% Revisions:      14/03/2013, small modifications and script description.
%
% Copyright (C) 2012 Tobias Siegfried
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details.
%% 7.1 Prepare
%profile on
clear all, clc
%cd('~/Dropbox/iMoMoPHASE1/model/iMoMoToolbox/enkf/models/rrModel/')
dataPath = '~/Dropbox/hydrosolutionsLocal/Projects_Active/SDC/iMoMoPHASE1/model/iMoMoToolboxTobi/app/themi/resources/trainingData/';
savePath = '~/Dropbox/hydrosolutionsLocal/Projects_Active/SDC/iMoMoPHASE1/model/iMoMoToolboxTobi/app/themi/resources/samples/'
saveName = 'themiSamples'
prepareLTTInputData('trainingData', dataPath, savePath, saveName)
%cd('~/Dropbox/iMoMoPHASE1/model/iMoMoToolbox/enkf/')
% 7.2 Run Assimilation (Training)
clc %clear all
tic
cd(savePath)
load themiSamples
cd('~/Dropbox/hydrosolutionsLocal/Projects_Active/SDC/iMoMoPHASE1/model/iMoMoToolboxTobi/src/enkf')

paths.home = '~/Dropbox/hydrosolutionsLocal/Projects_Active/SDC/iMoMoPHASE1/model/iMoMoToolboxTobi/src/enkf/';
paths.dataAssim = '~/Dropbox/hydrosolutionsLocal/Projects_Active/SDC/iMoMoPHASE1/model/data/processed/biascorrected/';
paths.prm = '~/Dropbox/hydrosolutionsLocal/Projects_Active/SDC/iMoMoPHASE1/model/iMoMoToolboxTobi/app/themi/prm/';
paths.forecasts = '~/Dropbox/hydrosolutionsLocal/Projects_Active/SDC/iMoMoPHASE1/model/data/processed/forecast';

[prm, x, x_true, E, stats, EStore] = fmain('../../app/themi/prm/prm-Themi-training.txt',[],E, [], paths);

cd('~/Dropbox/hydrosolutionsLocal/Projects_Active/SDC/iMoMoPHASE1/model/iMoMoToolboxTobi/app/themi/resources/trainingResults/');
clear x_true % not really the true state - address in new version!
save('themiTrained')
toc
% save training E
save([prm.customprm.locS0G0 prm.customprm.runType '/ETrain'],'E');
% at the same time

%profile off
%profile viewer
%% visualize
catch2Show = 2; 

% 7.3 Analyze Training (no snow)
a1 = squeeze(EStore(5 + (catch2Show-1)*8,:,:)); % Qcalc
a2 = squeeze(EStore(6 + (catch2Show-1)*8,:,:)); % Scalc
a3 = squeeze(EStore(7 + (catch2Show-1)*8,:,:)); % Gcalc
SMax = squeeze(EStore(8 + (catch2Show-1)*8,:,:)); % Smax

visStep = 50;

figure(15)
subplot(4,1,1)
boxplot(a1(:,1:visStep:end))
hold on
if ~prm.realWorld
    plot(prm.trueX(5 + (catch2Show-1)*8,1:visStep:end),'g')
    title(['true value: ' num2str(prm.trueX(5 + (catch2Show-1)*8)) ' (Rainfall retention parameter)']) 
else
    title('Rainfall retention parameter') 
end
hold off
xlabel('assim. step'),ylabel('/alpha_1')

grid minor
subplot(4,1,2)
boxplot(a2(:,1:visStep:end))
hold on
if ~prm.realWorld
    plot(prm.trueX(6 + (catch2Show-1)*8,1:visStep:end),'g')
    title(['true value: ' num2str(prm.trueX(6 + (catch2Show-1)*8)) ' (Evapotranspiration efficiency)']) 
else
    title('Evapotranspiration efficiency') 
end
hold off
xlabel('assim. step'),ylabel('/alpha_2')
grid minor
subplot(4,1,3)
boxplot(a3(:,1:visStep:end))
hold on
if ~prm.realWorld
    plot(prm.trueX(7 + (catch2Show-1)*8,1:visStep:end),'g')
    title(['true value: ' num2str(prm.trueX(7 + (catch2Show-1)*8)) ' (Groundwater storage discharge parameter)']) 
end
hold off
xlabel('assim. step'),ylabel('d')

grid minor
subplot(4,1,4)
boxplot(SMax(:,1:visStep:end))
hold on
if ~prm.realWorld
    plot(prm.trueX(8 + (catch2Show-1)*8,1:visStep:end),'g')
    title(['true value: ' num2str(prm.trueX(8 + (catch2Show-1)*8))]) 
end
hold off
xlabel('assim. step'),ylabel('Smax')

grid minor
%% VALIDATION RUN: Plot a run with the estimated parameters
% 17.12.2012, propagate whole ensemble for visualization of uncertainty!
% For the moment catchment no. 1
%cd('~/Dropbox/iMoMoPHASE1/model/iMoMoToolbox/enkf/models/rrModel/trainingData/themi')
%load('trainingData','themiSUB');

clc

catch2Show = 1;

prmV = prm;

lengthValPeriod = 2 * 365;

prmV.custom.serialDateStart = prm.current.serialDate + 1;
prmV.tValidation = (prmV.current.serialDate + 1: prmV.current.serialDate + 1 + lengthValPeriod)';

prmV.current.step = 1;
prmV.customprm.S0 = x(2:8:end)';
prmV.customprm.G0 = x(3:8:end)';
%prmV.customprm.SN0 = x(4);
prmV.customprm.limCutoff = 10^-6;
prmV.ensembleN = 100;
prmV.nstep = length(prmV.tValidation);
prmV.customprm.TmeanAvail = 1; % mean temperature available...
true_field = 0; 

idPStates = repmat((1:4)',prm.customprm.nC,1);
addID = repmat(1:prm.customprm.nC,4,1); addID = addID(:) - 1;
idPStates = idPStates + addID * prm.n;
        
idPar = repmat((5:8)',prm.customprm.nC,1);
addID = repmat(1:prm.customprm.nC,4,1); addID = addID(:) - 1;
idPar = idPar + addID * prm.n;

xO = E(idPStates,:);
p = E(idPar,:);
xFull = zeros(4 * prm.customprm.nC,prmV.nstep,prmV.ensembleN);

if ~prm.realWorld
    xOTrue = x(idPStates,1);
    pTrue = prm.trueX(idPar,1);
    xFullTrue = zeros(4 * prm.customprm.nC,prmV.nstep,1);
    prmV.trueX = [prm.trueX [xFullTrue; repmat(pTrue,1,size(xFull,2))]];
end

% start up phase
s = 1;
tFromBeg = prm.current.serialDate - prm.customprm.serialDateStart + 1;
tic
prmV.current.doy = datenum2doy(datevec(prmV.tValidation(prmV.current.step)));

xN = budykoStepFSVec(prmV,xO,p,true_field);

xFull(:,prmV.current.step,1:prmV.ensembleN) = xN; xO = xN;

if ~prm.realWorld
    xNTrue = budykoStepFSVec(prmV,xOTrue,pTrue,1);
    xFullTrue(:,prmV.current.step,1) = xNTrue; xOTrue = xNTrue;
    prmV.trueX(idPStates,tFromBeg) = xNTrue;
end
    
prmV.storage.SG(1:2:end,s,1:prmV.ensembleN) = xN(2 : 4 : end,:);
prmV.storage.SG(2:2:end,s,1:prmV.ensembleN) = xN(3 : 4 : end,:);
% move time
prmV.current.step =  prmV.current.step + 1;
prmV.current.serialDate =  prmV.current.serialDate + 1;
tFromBeg = tFromBeg + 1;

% remainder
for s = prmV.current.step : length(prmV.tValidation)
    
    if ~mod(s,10), s, end
    prmV.current.doy = datenum2doy(datevec(prmV.tValidation(prmV.current.step)));
    
    xN = budykoStepFSVec(prmV,xO,p,true_field);
    
    xFull(:,prmV.current.step,1:prmV.ensembleN) = xN; xO = xN;
    
    if ~prm.realWorld
        xNTrue = budykoStepFSVec(prmV,xOTrue,pTrue,1);
        xFullTrue(:,prmV.current.step,1) = xNTrue; xOTrue = xNTrue;
        prmV.trueX(idPStates,tFromBeg) = xNTrue;
    end
        
    prmV.storage.SG(1:2:end,prmV.current.step,1:prmV.ensembleN) = xN(2:4:end,:); % S
    prmV.storage.SG(2:2:end,prmV.current.step,1:prmV.ensembleN) = xN(3:4:end,:); % G
    % move time
    prmV.current.step =  prmV.current.step + 1;
    prmV.current.serialDate =  prmV.current.serialDate + 1;
    tFromBeg = tFromBeg + 1;
end
toc

%if prm.realWorld % write ETact observed in corresponding field
%    xFullTrue(4:4:end,:)  = ET(prm.current.step+1:prm.current.step+1+lengthValPeriod,:);
%end




%%
figure(1) % here, we only plot the the runoffs of the individual catchments - 
% x-checking the routing
for i = 1 : 8
    qConvFac = 1/(1000*24*3600) * prm.shp(prm.customprm.cID(i)).Area_m2;
    subplot(4,2,i)
    plot(prmV.tValidation,mean(qConvFac * xFull(1 + (i - 1) * 4,:,:),3),'r') 
    ylabel('[m^{3}/s]')
    datetick('x',10)
    title(strcat('Runoff Catchment ', num2str(i)))
    grid on
end

%%
figure(2)
hold on
cc=hsv(8);
for i = 1 : 8
    qConvFac = 1/(1000*24*3600) * prm.shp(prm.customprm.cID(i)).Area_m2;
    plot(prmV.tValidation,mean(qConvFac * xFull(1 + (i - 1) * 4,:,:),3),'color',cc(i,:))
    ylabel('[m^{3}/s]')
    datetick('x',10)
end
title('Runoff Themi Subcatchments')
grid minor
hold off
legend('C1','C2','C3','C4','C5','C6','C7','C8',2)

%%
figure(catch2Show+100)
qConvFac = 1/(1000*24*3600) * prm.shp(prm.customprm.cID(catch2Show)).Area_m2;

% runoff
subplot(4,1,1)
hold on
pol = [[prmV.tValidation;flipud(prmV.tValidation)] ...
    [min(qConvFac * xFull(1 + (catch2Show - 1) * 4,:,:),[],3)'; ...
    flipud(max(qConvFac * xFull(1 + (catch2Show - 1) * 4,:,:),[],3)')]];
fill(pol(:,1),pol(:,2),[0.9 .9 .9])
plot(prmV.tValidation,mean(qConvFac * xFull(1 + (catch2Show - 1) * 4,:,:),3),'r')
if ~prm.realWorld
    plot(prmV.tValidation,qConvFac * xFullTrue(1 + (catch2Show - 1) * 4,:),'b+')
end
hold off
legend('Ensemble spread','Q(modelled)','Q(real)',2)
title(strcat('Runoff Catchment ', num2str(catch2Show)))
ylabel('[m^{3}/s]')
datetick('x',10)
xlim([min(prmV.tValidation) max(prmV.tValidation)])
grid minor

% soil moisture
subplot(4,1,2)
hold on

pol = [[prmV.tValidation;flipud(prmV.tValidation)] ...
    [min(xFull(2 + (catch2Show - 1) * 4,:,:),[],3)'; ...
    flipud(max(xFull(2 + (catch2Show - 1) * 4,:,:),[],3)')]];

fill(pol(:,1),pol(:,2),[0.9 .9 .9])

plot(prmV.tValidation,mean(xFull(2 + (catch2Show - 1) * 4,:,:),3),'r')
if ~prm.realWorld
    plot(prmV.tValidation,xFullTrue(2 + (catch2Show - 1) * 4,:),'b+')
end
hold off
legend('Ensemble spread','S(modelled)','S(real)',2)
title('Soil Moisture')
ylabel('[mm]')
datetick('x',10)
xlim([min(prmV.tValidation) max(prmV.tValidation)])
grid minor

% groundwater
subplot(4,1,3)
hold on

pol = [[prmV.tValidation;flipud(prmV.tValidation)] ...
    [min(xFull(3 + (catch2Show - 1) * 4,:,:),[],3)'; ...
    flipud(max(xFull(3 + (catch2Show - 1) * 4,:,:),[],3)')]];

fill(pol(:,1),pol(:,2),[0.9 .9 .9])

plot(prmV.tValidation,mean(xFull(3 + (catch2Show - 1) * 4,:,:),3),'r')
if ~prm.realWorld
    plot(prmV.tValidation,xFullTrue(3 + (catch2Show - 1) * 4,:),'b+')
end
hold off
legend('Ensemble spread','G(modelled)','G(real)',2)
title('Groundwater')
ylabel('[mm]')
datetick('x',10)
xlim([min(prmV.tValidation) max(prmV.tValidation)])
grid minor

% actual evaporation
subplot(4,1,4)
hold on
pol = [[prmV.tValidation;flipud(prmV.tValidation)] ...
    [min(xFull(4 + (catch2Show - 1) * 4,:,:),[],3)'; ...
    flipud(max(xFull(4 + (catch2Show - 1) * 4,:,:),[],3)')]];
fill(pol(:,1),pol(:,2),[0.9 .9 .9])
plot(prmV.tValidation,mean(xFull(4 + (catch2Show - 1) * 4,:,:),3),'r')
if ~prm.realWorld
    plot(prmV.tValidation,xFullTrue(4 + (catch2Show - 1) * 4,:),'b+')
end
hold off

legend('Ensemble spread','Etact(modelled)','ETact(real)',2)
title('Evapotranspiration (actual)')
ylabel('[mm/day]')
datetick('x',10)
xlim([min(prmV.tValidation) max(prmV.tValidation)])
grid minor