function prepareLTTInputData(fName, dataPath, savePath, saveName)
%% function prepareLTTInputData(fName, datapath, savePath, saveName)
%
% This function prepares the long-term training data and stores it in the 
% appropriate location. 
%
% @param fName      - Loads basin-specific training data
% @param dataPath   - Path where basin-specific training data is stored.
% @param savePath   - Output path
% @param saveName   - Output name
%
% Usage:
%                 prepareLTTInputData(fName, datapath, savePath, saveName)
%
% File:           prepareLTTInputData.m
%
% Created:        06/12/2012
%
% Last modified:  30/04/2013
%
% Author:         Tobias Siegfried (hydrosolutions ltd.)           
%
% Purpose:        Model training preparation.
%
% Description:    
%
% Data preparation function. It should be noted that this
% function is not yet adapted to the case where snow is present in a basin.
%
% Revisions:      30/04/2013, 
%                 14/03/2013, function help tweaked and expanded
%
% Copyright (C) 2012 Tobias Siegfried
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details.
%% Go to data location
cd(dataPath)
load(fName)
%% Ensemble preparation
% random number generator
nG = rng;

n = 8; % n is state parameter space (this is also in the prm-Themi.txt file)
n_sample = 2000; %number of samples for the ensemble (ditto)
effPor = 0.1;

c = [1 2 3 4 5 6 7 8]; % id of catchments to b used in assimilation run (right now, we have 8 catchments).
cID = c;

% states (ensembles)
qS = 10 + 2 * randn(length(c),n_sample); % 'assumed' mean runoff - specific runoff [mm]
SS =  emprand(mean(SMTrainSub * effPor,2),length(c),n_sample); % accounting for effective porosity
SG = 10 + 2 * randn(length(c),n_sample); % assumed groundwater ensemble - UNITS
SSN = zeros(length(c),n_sample); % UNITS (absolutely no snow here in Themi)
SET = emprand(mean(ETactTrainSub,2),length(c),n_sample);

% parameters
alpha1 = betarnd(2,2,[length(c) n_sample]);             % 3 alphas
alpha2 = betarnd(2,2,[length(c) n_sample]);
d = betarnd(2,2,[length(c) n_sample]);

p4 = max(SS(:)) + 50 * randn(length(c),n_sample);       % Smax
p5 = 0.1 * randn(length(c),n_sample);                   % rain-snow transition temperature. Note: the spread around the zero mean is entirely guessed here!
p6 = 4.5 + 0.6 * randn(length(c),n_sample);             % DD melting factor in mm H2O / day / degC

% full ensemble (with snow)
%E = [qS;SS;SG;SSN;SET;alpha1;alpha2;d;p4;p5;p6]; valid for a single catchment
E = cat(3,qS', SS', SG', SSN', SET', alpha1', alpha2', d', p4', p5', p6');

% full ensemble (without snow)
% E = [qS;SS;SG;SET;alphaS;p4]; % valid for a single catchment
E = cat(3,qS', SS', SG', SET', alpha1', alpha2', d', p4');

E = shiftdim(E,2);

% x_true (snow)
% x_true = [rand(5,1); 0.2714; 0.5156; 0.0188; 700.0950;0;.5];

% x_true (no snow)
% x_true = [rand(4,length(c)); repmat([0.2714; 0.5156; 0.0188; 32],1,length(c))];
x_true = [rand(4,length(c)); rand(4, length(c))]; x_true(end,:) = x_true(end,:) * 30;


% Training data
ET = ETactTrainSub(:,c);
P = PTrainSub(:,c);
Tmax = TmaxTrainSub(:,c);
Tmin = TminTrainSub(:,c);
Tmean = TmTrainSub(:,c);

SM = SMTrainSub * effPor; % accounting for effective porosity

Lat = nanmean(themiSUB(1).Y); % only first catchment

% here, assume that there is no snow (just to constrain the model better)!
SN = SMTrainSub * 0;

% other parameters
model = 'rrModel';
t = 0;
time = tTraining;
seed = 1;

% inital catchment storage is 0
S0 = zeros(1,length(c));
G0 = S0;
%% saving the samples (this should go into the ~/samples/ path)
cd(savePath)
save(saveName,'n','n_sample','x_true','E','model','seed','t','time','P','ET','Tmin','Tmax','Tmean','SM','Lat','SN','cID');

% this is a bit of hardwiring - adress in update
cd('../restart/training/')
trainingStart = datevec(now);

fileID = fopen('simTime.txt','w');
fprintf(fileID,'%e %e %e %e %e %e',trainingStart);
fclose(fileID);
save('S0G0','S0','G0')
end

