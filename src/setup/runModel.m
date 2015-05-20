 %% Update and run model
% Script to process the raw data from the internet.
%
% File:          runModel.m
%
% Created:        05/02/2015
%
% Last modified:  05/02/2015
%
% Author:         Sebastian Stoll (hydrosolutions ltd.)
%
% Purpose:        Script that assimilates the states and run the model.
%
% Description:   Script that assimilates the states and run the model.
%
%
% Copyright (C) 2015 hydrosolutions
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details. 
 
 
 
try 
   tic
    %% SETTINGS
    % LOAD SETUP FILE
    clear all
    
    
    load('setup.mat')
    
    sPath = setup.Path(1:end-4);
    
    % COLLECTION OF PATHS (FOR ASSIMILATION)
    paths.home = strcat(sPath,'/src/enkf/');
    paths.main = strcat(setup.mPath,'/');
    paths.data = strcat(setup.mPath,'/data/');
    
    %% 1. RUN MODEL    
    % Run Model
    cd(paths.home)
    assimilation(paths);
    
    
    %% 2. PRODUCE DATABASE DATA
    cd(strcat(paths.data,'processed/sub/'))
        
    % METEO DATA
    % Precipitation
    load('SUB_FP.mat');
    Pmax=max(FPsub{end}.FP,[],3);
    Pmin=min(FPsub{end}.FP,[],3);
    Pmean=mean(FPsub{end}.FP,3);
    Ptime=FPsub{end}.timeF;
    dateMet=cellstr(strcat(datestr(Ptime,'yyyy'),'-',datestr(Ptime,'mm'),'-',datestr(Ptime,'dd'),' 12:00:00'));
    % Mean, max, min
    P(:,:,1)=Pmean;
    P(:,:,2)=Pmax;
    P(:,:,3)=Pmin;

    % Temperature
    load('SUB_FT.mat');
    Tmax=max(FTsub{end}.FmaxT,[],3);
    Tmin=min(FTsub{end}.FminT,[],3);
    Tmean=mean(FTsub{end}.FT,3);
    % Mean, max, min
    T(:,:,1)=Tmean;
    T(:,:,2)=Tmax;
    T(:,:,3)=Tmin;
    
    %MODEL SIMULATIONS
    
   
    % Q
    cd(strcat(paths.main,'results/forecast/Q'))
    fn = dir;
    fn = {fn(end).name}';
    fn = cell2mat(fn);
    load(fn);

    % S
    cd(strcat(paths.main,'results/forecast/S'))
    fn = dir;
    fn = {fn(end).name}';
    fn = cell2mat(fn);
    load(fn);

    % G
    cd(strcat(paths.main,'results/forecast/G'))
    fn = dir;
    fn = {fn(end).name}';
    fn = cell2mat(fn);
    load(fn);
    
    % ET
    cd(strcat(paths.main,'results/forecast/ET'))
    fn = dir;
    fn = {fn(end).name}';
    fn = cell2mat(fn);
    load(fn);

    for i=1:size(Pmean,2) % subcatchments
    % Discharge
        % model
        for j=1:length(FQ.catch{i}.Q)
            maxMQ(:,j)=max(FQ.catch{i}.Q{j},[],2);
            minMQ(:,j)=min(FQ.catch{i}.Q{j},[],2);
            meanMQ(:,j)=mean(FQ.catch{i}.Q{j},2);
        end
        % ensemble
        Q(:,i,2)=max(maxMQ,[],2);
        Q(:,i,3)=min(minMQ,[],2);
        Q(:,i,1)=mean(meanMQ,2);

    % Soil moisture
        % model
        for j=1:length(FS.catch{i}.Q)
            maxMS(:,j)=max(FS.catch{i}.Q{j},[],2);
            minMS(:,j)=min(FS.catch{i}.Q{j},[],2);
            meanMS(:,j)=mean(FS.catch{i}.Q{j},2);
        end
        % ensemble
        S(:,i,2)=max(maxMS,[],2);
        S(:,i,3)=min(minMS,[],2);
        S(:,i,1)=mean(meanMS,2);

    % Groundwater
        % model
        for j=1:length(FG.catch{i}.Q)
            maxMG(:,j)=max(FG.catch{i}.Q{j},[],2);
            minMG(:,j)=min(FG.catch{i}.Q{j},[],2);
            meanMG(:,j)=mean(FG.catch{i}.Q{j},2);
        end
        % ensemble
        G(:,i,2)=max(maxMG,[],2);
        G(:,i,3)=min(minMG,[],2);
        G(:,i,1)=mean(meanMG,2);

    % Evapotranspiration
        % model
        for j=1:length(FET.catch{i}.Q)
            maxMG(:,j)=max(FET.catch{i}.Q{j},[],2);
            minMG(:,j)=min(FET.catch{i}.Q{j},[],2);
            meanMG(:,j)=mean(FET.catch{i}.Q{j},2);
        end
        % ensemble
        ET(:,i,2)=max(maxMG,[],2);
        ET(:,i,3)=min(minMG,[],2);
        ET(:,i,1)=mean(meanMG,2);    

    end


    %Simulation date
    dateSim=cellstr(strcat(datestr(FQ.timeF,'yyyy'),'-',datestr(FQ.timeF,'mm'),'-',datestr(FQ.timeF,'dd'),' 12:00:00'));

    % save matrix
    cd(strcat(paths.data,'processed/db/'))
    today=num2str(datenum(date));
    save(today,'P', 'T', 'Q', 'S', 'G', 'ET', 'dateMet', 'dateSim');
    clear all
    toc
    
catch
    disp('Error occured, next try in 3 hours')
end
    
    
    
    
    
    
    
    
    
    
    
    
    