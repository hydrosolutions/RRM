function out = assimilation(paths)
%% function out = onlineAssiModel(warmup,paths,choiceE)
%
% This function performs the daily assimilation step.
%
% % @param paths      - loads all paths from an external path file (see pathsFile.m)
% % @out              - returns 1 in case of success.
%
% Usage:
%                 out = assimilation(paths)
%
% File:           assimilation.m
%
% Created:        05/01/2013
%
% Last modified:  02/02/2015
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

cd(strcat(paths.main,'/prm/'))
nN = 35;
fid = fopen('prm-Global.txt');
parmN = textscan(fid,'%s',nN,'CommentStyle','%');
parmS = textscan(fid,'%s','CommentStyle','%');
fclose(fid);
cellfun(@evalc, parmN{1},'UniformOutput',0);
cellfun(@evalc, parmS{1},'UniformOutput',0);


%% READ CUSTOM PARAMETERS
nN = 19;
fid = fopen('prm-Custom.txt');
parmN = textscan(fid,'%s',nN,'CommentStyle','%');
parmS = textscan(fid,'%s','CommentStyle','%');
fclose(fid);
cellfun(@evalc, parmN{1},'UniformOutput',0);
cellfun(@evalc, parmS{1},'UniformOutput',0);
prm.customprm = customprm;


%% CHECK AVAILABILITY OF ASSIMILATION DATA & DATA PREPARATION
cd(strcat(paths.main,'/data/processed/sub/'))
load('SUB_P'); % FEWS 0.1 deg resolution, precipitation
load('SUB_pET'); % Based on NCEP CFS Temperature data, using Hargreaves, GDAS
load('SUB_FP'); % NCEP GFS (precipitation forecasts)
load('SUB_FpET'); % NCEP GFS (ET forecasts)

% load(strcat(paths.main,'/resources/geometry/rout.mat')); % routing information for subcatchments
currDate = datevec(FPsub{end}.timeF(1));


%% LOAD PREEXISTING DATA
cd(strcat(paths.main,'/prm/'));
[prm] = model_setprm(prm);
prm.path.E=strcat(paths.main,'resources/restart/E.mat');
prm.path.S0G0=strcat(paths.main,'resources/restart/S0G0.mat');

load(prm.path.E,'E');
Ein = reshape(E, [prm.n size(E,2) size(E,1)/prm.n]); % just get Ein (E0) in shape!

% general prep.
prm.customprm.nC = length(prm.customprm.cID);
for i = 1 : prm.customprm.nC
    prm.customprm.qConvFac(i) = 1/(1000*24*3600) * prm.shp(prm.customprm.cID(i)).Area; % m3/s
end

%% Assimilation
% Note that assimilation can only be carried out based on data from 'yesterday' and only when data from
% yesterday is acutally complete in the sense that GDAS (Temp and then ET) and FEWS (Precip.) data are
% available. Since FEWS is only becoming available late night (CET), we
% might utlimately just run the one day assimilation step around midnight CET.

%Settings
prm.ensembleN = prm.m;
prm.current.step = 1;
prm.nstep = 1;
prm.customprm.simStartY = currDate(1); % this corresponds to the first forecast day.
prm.customprm.simStartM = currDate(2);
prm.customprm.simStartD = currDate(3);
prm.customprm.time = timeP(end);
prm.customprm.P = Psub(end,:);
prm.customprm.ET = pETsub(end,:);

load(prm.path.S0G0,'SG');
prm.storage.S0G0 = SG;

%Check for Discharge data
Q = getObservation('Q', prm.customprm.time,prm.customprm.nC,paths); 

%Assign Values
prm.customprm.Q = Q;
prm.customprm.Tmean = NaN(size(prm.customprm.P));
prm.customprm.Tmin = NaN(size(prm.customprm.P));
prm.customprm.Tmax = NaN(size(prm.customprm.P));
prm.customprm.SM = NaN(size(prm.customprm.P));
prm.customprm.G = NaN(size(prm.customprm.P));


% ONE STEP ASSIMILATION (USING DATA FROM 'YESTERDAY')
cd(paths.home)
[prm, x, x_true, Eout, stats, EStore] = main(prm,[],Ein);

% PROCESS / STORE UPDATED STATES ETC.
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
    prm.customprm.ET = squeeze(FpETsub{end}.FpET(:,:,ensF));
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
sCatch = 3;
prm.forecasts.T = FPsub{end}.timeF;
prm.forecasts.P = squeeze(FPsub{end}.FP(:,sCatch,:));
prm.forecasts.ET = squeeze(FpETsub{end}.FpET(:,sCatch,:));

%[fH1,fH2,fH3,fH4,fH4] = plotEnsemble(prm,xFullF,sCatch,figNum,holdSP);

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


cd(strcat(paths.main,'results/forecast/'))
cd('Q'), save([num2str(prm.current.SerialDateVec(1)) '_FQ'],'FQ'),
cd('..')
cd('S'), save([num2str(prm.current.SerialDateVec(1)) '_FS'],'FS')
cd('..')
cd('G'), save([num2str(prm.current.SerialDateVec(1)) '_FG'],'FG'),
cd('..')
cd('ET'), save([num2str(prm.current.SerialDateVec(1)) '_FET'],'FET'),
cd('..')

% Save nowcasts
cd(strcat(paths.main,'results/nowcast/'))
save([num2str(prm.customprm.time(1)) '_E'],'E')


%% FEEDBACK
out = 1;

end



