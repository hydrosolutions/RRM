function out = onlineAssiModel(warmup,paths,choiceE)
%% function out = onlineAssiModel(warmup,paths,choiceE)
%
% This function performs the daily assimilation step.
%
% @param prmfname   - Loads basin-specific training data (i.e. for Themi,'/prm/prm-Themi.txt')
% @param warmup     - perform a warmup run with data stored data from GDAS.
% @param paths      - loads all paths from an external path file (see pathsFile.m)
% @param choiceE    - selection of E matrix (either 'training' or 'online')
% @out              - returns 1 in case of success.
%
% Usage:
%                 out = onlineAssiModel(warmup,paths,choiceE)
%
% File:           onlineAssimilationModel.m
%
% Created:        05/01/2013
%
% Last modified:  27/11/2014
%
% Author:         Tobias Siegfried (hydrosolutions ltd.)
%
% Purpose:        Online Assimilation
%
% Description:    Wrapper for online data assimilation and flow forecasting.
% For the iMoMo application, prmfname = 'prm/prm-Themi.txt'. For other
% applications, change accordingly and make sure that models and prm files
% are stored in the appropriate locations. 'warmup' gives the option of
% 'bridging the training and the online assimilation periods with an
% intermediate warmup.
%
% Data preparation function. It should be noted that this
% function is not yet adapted to the case where snow is present in a basin.
%
% Revisions:      27/11/2014, simplified parameter reading and file
%                 handling. Simplification of path structure.
%                 08/05/2013, adding the option of choosing E for starters
%                 07/05/2013, revising issues
%                 14/03/2013, function help tweaked and expanded
%                 11 / 02 / 2013, addressing the S0 & G0 issue.
%                 04/02/2013, first 'real' online forecasts. S0 and G0
%                 issue as well as pre-wetting still outstanding.
%                 28/01/2013, work continues towards doing the online
%                 assimilation. One consideration to make is how to bring in best the
%                 'learned' parameter values from the training runs. sort that out first.
%
% Copyright (C) 2012 hydrosolutions ltd., Lindenbachstrasse 11, CH-8006
% Zurich, Switzerland
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details.
%% ASSUME ERROR
out = 0;
%% READ GLOBAL PARAMETERS
keyboard

cd(paths.prm)
nN = 35;
fid = fopen(paths.fileN.globalPar);
parmN = textscan(fid,'%s',nN,'CommentStyle','%');
parmS = textscan(fid,'%s','CommentStyle','%');
fclose(fid);
cellfun(@evalc, parmN{1},'UniformOutput',0);
cellfun(@evalc, parmS{1},'UniformOutput',0);
%% READ CUSTOM PARAMETERS
nN = 19;
fid = fopen(paths.fileN.modelPar);
parmN = textscan(fid,'%s',nN,'CommentStyle','%');
parmS = textscan(fid,'%s','CommentStyle','%');
fclose(fid);
cellfun(@evalc, parmN{1},'UniformOutput',0);
cellfun(@evalc, parmS{1},'UniformOutput',0);
prm.customprm = customprm;
%% CHECK AVAILABILITY OF ASSIMILATION DATA & DATA PREPARATION
cd(paths.data)
load('SUB_P'); % FEWS 0.1 deg resolution, precipitation
load('SUB_pET'); % Based on NCEP CFS Temperature data, using Hargreaves, GDAS
cd('..') % Right now, there is no bias correction for NCEP GFS (can become available once station data is coming online)
load('SUB_FP'); % NCEP GFS (precipitation forecasts)
load('SUB_FpET'); % NCEP GFS (ET forecasts)
load('rout'); % routing information for subcatchments
currDate = datevec(FPsub{end}.timeF(1));
%% LOAD PREEXISTING DATA
cd(paths.prm)
[prm] = model_setprm(prm);

if strcmp(choiceE,'online')
    load([prm.customprm.locS0G0 prm.customprm.runType '/E'],'E');
elseif strcmp(choiceE,'training')
    load([prm.customprm.locS0G0 'training/ETrain'],'E');
else
    error('No suitable E found. Please either choose ''online'' or ''training''.')
end

Ein = reshape(E, [prm.n size(E,2) size(E,1)/prm.n]); % just get Ein (E0) in shape!
%% MODEL WARMUP (WETTING)
% Note: Warmup uses last training SG and E for starters. Warmup can only be
% done until t-2 (if today date is t).

% general prep.
prm.customprm.nC = length(prm.customprm.cID);
for i = 1 : prm.customprm.nC
    prm.customprm.qConvFac(i) = 1/(1000*24*3600) * prm.shp(prm.customprm.cID(i)).Area; % m3/s
end

if warmup  
    % 0. copy training files (S0G0)
    cd(prm.customprm.locS0G0)
    cd('training')
    copyfile('S0G0.mat',[prm.customprm.locS0G0 prm.customprm.runType '/']);
    copyfile('simTime.txt',[prm.customprm.locS0G0 prm.customprm.runType '/']);
    cd(paths.home)
    % tag ensemble mean GFS{end-1} P data to FEWS so as to bridge the one day gap.
    timeP = [timeP FPsub{end-1}.timeF(1)];
    tagP = mean(FPsub{end-1}.FP(1,:,:),3);
    Psub = [Psub; tagP];
    % 1. just make sure that times are right and that the overlap between available P/ET warmup data
    [overlapT,iP,iT] = intersect(timeP,timeT);
    % Above Note: a bold way to overcome the one day gap due to lacking FEWS data is to
    % simply fill in with GFS data from yesterday and the first day
    % forecast data (i.e. from yesterday) for P.
    warmupP = Psub(iP,:);
    warmupET = pETsub(iT,:);
    prm.ensembleN = prm.m;
    % 2. write corresponding values into places.
    prm.customprm.time = timeP(iP)';
    prm.customprm.serialDateStart = prm.customprm.time(1);
    prm.current.step = 1;
    prm.current.SerialDateVec = datevec(prm.customprm.time(1));
    prm.current.serialDate = prm.customprm.serialDateStart;
    prm.customprm.P = warmupP;
    prm.customprm.ET = warmupET;
    % 2.1. just make sure that no bs remains around
    prm.customprm.Tmin = NaN(size(warmupP));
    prm.customprm.Tmax = NaN(size(warmupP));
    prm.customprm.Tmean = NaN(size(warmupP));
    prm.customprm.SM = NaN(size(warmupP));
    prm.customprm.Q = NaN(size(warmupP));
    prm.customprm.G = NaN(size(warmupP));
    prm.nstep = length(iP);
    prm.storage.SG = zeros(prm.customprm.nC * 2, prm.nstep, prm.ensembleN);
    % 3. get S0 and G0
    cd(prm.customprm.locS0G0)
    cd('online')
    load([prm.customprm.locS0G0 prm.customprm.runType '/S0G0.mat'], 'SG')
    prm.storage.S0G0 = SG;
    % 4. states & parameters
    idPStates = repmat((1:4)',prm.customprm.nC,1);
    addID = repmat(1:prm.customprm.nC,4,1); addID = addID(:) - 1;
    idPStates = idPStates + addID * prm.n;
    idPar = repmat((5:8)',prm.customprm.nC,1);
    addID = repmat(1:prm.customprm.nC,4,1); addID = addID(:) - 1;
    idPar = idPar + addID * prm.n;
    xO = E(idPStates,:);
    p = E(idPar,:);
    xFull = zeros(4 * prm.customprm.nC,prm.nstep,prm.ensembleN);
    for t = 1 : prm.nstep
        % do the steps
        prm.current.doy = datenum2doy(datevec(prm.current.serialDate));
        xN = budykoStepFSVec(prm,xO,p,0); % actual water balance
        xFull(:,prm.current.step,1:prm.ensembleN) = xN; xO = xN;
        prm.storage.SG(1:2:end,prm.current.step,1:prm.ensembleN) = xN(2:prm.p:end,:); % S
        prm.storage.SG(2:2:end,prm.current.step,1:prm.ensembleN) = xN(3:prm.p:end,:); % G
        % move time
        prm.current.step =  prm.current.step + 1;
        prm.current.serialDate =  prm.current.serialDate + 1;
        prm.current.serialDateVec = datevec(prm.current.serialDate);
    end
    SG = squeeze(prm.storage.SG(:,end,:));
    save([prm.customprm.locS0G0 prm.customprm.runType '/S0G0'],'SG','-append');
    % NOTE: no update of E necessary since we do not run any assimilation here!
end % if warmup
%% NEW OBSERVATIONS? - If yes, do assimilation, if no, simply run ensemble forecasts
% Note that assimilation can only be carried out based on data from 'yesterday' and only when data from
% yesterday is acutally complete in the sense that GDAS (Temp and then ET) and FEWS (Precip.) data are
% available. Since FEWS is only becoming available late night (CET), we
% might utlimately just run the one day assimilation step around midnight CET.
prm.ensembleN = prm.m;
newObservations = 0;
%if (timeP(end) + timeT(end) + 2)/2 == floor(now) % can be solved more elegantly later inside the timer function.
if ~newObservations % for the moment once all works, use line-above statement
    newObservations = 1; % there are new observations, at least from the satellite (note: only after 6pm CET available).
    % -> iMoMo DATABASE SCREENING
    prm.current.step = 1;
    prm.nstep = 1;
    prm.customprm.simStartY = currDate(1); % this corresponds to the first forecast day.
    prm.customprm.simStartM = currDate(2);
    prm.customprm.simStartD = currDate(3);
    prm.customprm.time = timeP(end);
    prm.customprm.P = Psub(end,:);
    load([prm.customprm.locS0G0 prm.customprm.runType '/S0G0'],'SG');
    prm.storage.S0G0 = SG;
    % T or ET
    if prm.customprm.calculateET
        % NA for Themi, fill in if required for other basins.
    else
        prm.customprm.ET = pETsub(end,:);
        prm.customprm.Tmean = NaN(size(prm.customprm.P));
        prm.customprm.Tmin = NaN(size(prm.customprm.P));
        prm.customprm.Tmax = NaN(size(prm.customprm.P));
        prm.customprm.SM = prm.customprm.Tmean;
        prm.customprm.Q = prm.customprm.Tmean;
        prm.customprm.G = prm.customprm.Tmean;
    end
else
    disp('Assimilation not necessary since no new data is available!')
    % to be worked on! TS 28/11/2014 - ???
end
%% ONE STEP ASSIMILATION (USING DATA FROM 'YESTERDAY')
if newObservations
    cd(paths.home)
    [prm, x, x_true, Eout, stats, EStore] = main(prm,[],Ein);
end

%keyboard
%% PROCESS / STORE UPDATED STATES ETC.
% if calculation, successful write back current values to 'prm/prm-custom-themi-online.txt'
[E] = write_customprm(paths,prm,Eout);
%% FLOW FORECAST USING UPDATED STATES
% this is now 5 days x n subcatchments x 21 ensemble runs
% NOTE: here, we can either loop over the individual forecast ensembles or,
% maybe even smarter, tag all the ensembles together and do the flow
% accounting correspondingly...

% just setting up basics
prm.nstep = prm.customprm.fclength; % forecast horizon
prm.ensembleNF = size(FPsub{end}.FP,3); % this is the number of ensemble forecasts
true_field = 0;
prm.storage.SG = zeros(prm.customprm.nC * 2, prm.nstep, prm.ensembleN);
prm.storage.SGF = repmat(prm.storage.SG,[1 1 1 prm.ensembleNF]);

% states & parameters
idPStates = repmat((1:4)',prm.customprm.nC,1);
addID = repmat(1:prm.customprm.nC,4,1); addID = addID(:) - 1;
idPStates = idPStates + addID * prm.n;
idPar = repmat((5:8)',prm.customprm.nC,1);
addID = repmat(1:prm.customprm.nC,4,1); addID = addID(:) - 1;
idPar = idPar + addID * prm.n;
xO = E(idPStates,:);
p = E(idPar,:);
xFull = zeros(4 * prm.customprm.nC,prm.nstep,prm.ensembleN);
xFullF = repmat(xFull, [1 1 1 prm.ensembleNF]);

for ensF = 1 : prm.ensembleNF % these are # ensemble forecasts
    % copy P / ET over to prm.customprm
    prm.customprm.P = squeeze(FPsub{end}.FP(:,:,ensF));
    prm.customprm.ET = squeeze(FpETsub{end}.pET(:,:,ensF));
    % time
    prm.current.step = 1; prm.current.SerialDateVec = FPsub{end}.timeF;
    prm.customprm.serialDateStart = FPsub{end}.timeF(1);
    prm.current.serialDate = prm.customprm.serialDateStart;
    % actual forecast run
    for tf = 1 : 5 % this is the GFS forecast horizon, i.e. 5 days.
        prm.current.doy = datenum2doy(datevec(FPsub{end}.timeF(tf)));
        xN = budykoStepFSVec(prm,xO,p,true_field);
        xFull(:,prm.current.step,1:prm.ensembleN) = xN; xO = xN;
        prm.storage.SG(1:2:end,prm.current.step,1:prm.ensembleN) = xN(2:4:end,:); % S
        prm.storage.SG(2:2:end,prm.current.step,1:prm.ensembleN) = xN(3:4:end,:); % G
        % move time
        prm.current.step =  prm.current.step + 1;
        prm.current.serialDate =  prm.current.serialDate + 1;
    end
    % store results correspondingly (now, in 4D matrices).
    prm.storage.SGF(:,:,:,ensF) = prm.storage.SG;
    xFullF(:,:,:,ensF) = xFull;
end

%% PROCESS (AND VISUALIZE) RESULTS AND STORE RELEVANT DATA
holdSP = 1; sCatch = 3; figNum = 1;

prm.forecasts.T = FPsub{end}.timeF;
prm.forecasts.P = squeeze(FPsub{end}.FP(:,sCatch,:));
prm.forecasts.ET = squeeze(FpETsub{end}.pET(:,sCatch,:));

% [fH1,fH2,fH3,fH4,fH4] = plotEnsemble(prm,xFullF,sCatch,figNum,holdSP);

% store data in /Dropbox/.../model/data/processed/forecast/...
for c = 1 : prm.customprm.nC
    qConvFac = 1;
    for ensF = 1 : prm.ensembleNF
        FQ.timeF = prm.current.SerialDateVec;
        FQ.catch{c}.Q{ensF} = qConvFac * squeeze(xFullF(1 + (c - 1) * 4,:,:,ensF));
        FQ.units = 'm3/s';
        FS.timeF = prm.current.SerialDateVec;
        FS.catch{c}.Q{ensF} = squeeze(xFullF(2 + (c - 1) * 4,:,:,ensF));
        FS.units = 'mm';
        FG.timeF = prm.current.SerialDateVec;
        FG.catch{c}.Q{ensF} = squeeze(xFullF(3 + (c - 1) * 4,:,:,ensF));
        FG.units = 'mm';
        FET.timeF = prm.current.SerialDateVec;
        FET.catch{c}.Q{ensF} = squeeze(xFullF(4 + (c - 1) * 4,:,:,ensF));
        FET.units = 'mm';
    end
end

cd(paths.forecasts)
cd('Q'), save([num2str(prm.current.SerialDateVec(1)) '_FQ'],'FQ'),
cd('..')
cd('S'), save([num2str(prm.current.SerialDateVec(1)) '_FS'],'FS')
cd('..')
cd('G'), save([num2str(prm.current.SerialDateVec(1)) '_FG'],'FG'),
cd('..')
cd('ET'), save([num2str(prm.current.SerialDateVec(1)) '_FET'],'FET'),
cd('..')
%% FEEDBACK
out = 1;
if out
    disp('============================================')
    disp('Online Assimilation Successfully Terminated!')
    disp('============================================')
end

keyboard

end



