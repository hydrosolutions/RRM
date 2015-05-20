function processRaw
% Script to process the raw data from the internet.
%
% File:           getRaw.m
%
% Created:        01/12/2014
%
% Last modified:  01/12/2014
%
% Author:         Sebastian Stoll (hydrosolutions ltd.)
%
% Purpose:        Process raw data and save them.
%
% Description:    Process raw data and save them.
%
%
% Copyright (C) 2014 hydrosolutions
%
% This file is part of iMoMo-Matlab. iMoMo-Matlab is a free software and
% licensed under the Free Software Foundation. See LICENSE for details.


try

%% SETUP
tic
clear all

warning off

%Today's date
today=datenum(date);

%load setup file
load('setup.mat')

%bounding box
bblon = [setup.bBox(1) setup.bBox(3)];
bblat = [setup.bBox(2) setup.bBox(4)];
ndays = 0; %number of days for download of GDAS/GFS (Maximum = 6)-1
pfad = setup.mPath;
region = setup.ndvi; %in future read this information from file
cd(pfad)


    %% 1.Load NDVI
    cd('./data/raw/ndvi')
    if strcmp(setup.ndvi,'N/A')==0
        %indentify and read latest NDVI image
        fn = dir;
        fn = {fn(end).name}';
        fn = cell2mat(fn);
        load(fn);

        %convert to RGB to NDVI
        NDVI=single(NDVI)/250;

        %resize
        [NDVI,R_NDVI] = resizem(NDVI,0.5,R_NDVI);

        %convert to Kc
        Kc=1.4571*(NDVI)-0.1725; %  according to: Kamble, B., Irmak, A., & Hubbard, K. (2013). Estimating Crop Coefficients Using Remote Sensing-Based Vegetation Index. Remote Sensing, 5(4), 1588?1602. doi:10.3390/rs5041588

        %set clouds and open water to 1
        Kc(Kc<=0)=1;

        %clear variables
        clearvars fn NDVI
    end
    
    %% 2. Adjust  GDAS T to elevation
    % load GDAS data
    cd('..')
    load(strcat(setup.T,'.mat'));
            
    % specify parameter
    param{2}.lapse=-0.5;
    
    %check if DEM is available
    cd('./dem')
    if exist(fullfile(cd, 'DEM_GDAS.mat'), 'file')==0
        % GET DEM FROM NASA WORLDWIND
        [DEM_GDAS,R_DEM_GDAS] = getDEM(R_T.LatitudeLimits, R_T.LongitudeLimits, 10);
        DEM_GDAS=double(DEM_GDAS);
        % Resize and sample
        [DEM_GDAS,R_DEM_GDAS] = resizem(DEM_GDAS,[size(T,1) size(T,2)],R_DEM_GDAS);
        % Save
        save 'DEM_GDAS.mat' DEM_GDAS R_DEM_GDAS
    else
        load('DEM_GDAS.mat')
    end
           
    %check if nearest neighbor for GDAS raster has already been calculated
    if exist(fullfile(cd, 'iEiC_GDAS.mat'), 'file')==0
            [iE_GDAS, iC_GDAS] = findNeighbor(1500,R_T,DEM_GDAS);
            save iEiC_GDAS.mat iE_GDAS iC_GDAS;
    else
        load('iEiC_GDAS.mat');
    end
                      
    %Adjust Temperatures
    [T] = adjMeteodata(T, param, 2, DEM_GDAS, iC_GDAS, iE_GDAS); % adjust to elevation
    [maxT] = adjMeteodata(maxT, param, 2, DEM_GDAS, iC_GDAS, iE_GDAS); % adjust to elevation
    [minT] = adjMeteodata(minT, param, 2, DEM_GDAS, iC_GDAS, iE_GDAS); % adjust to elevation
    
    cd('..');

    
    %% 3. Bias correct GDAS T
    
    %Load station data
    load WMO.mat;
    statData.WMO = WMO;
    statData.iMoMo = WMO;
    
    % Specify 'believe' weigths that are an indication on how much data can be trusted (NOTE: sum(weigths) = 1 required).
    statData.w.WMO = 1;
    statData.w.iMoMo = 0;
    
    % Specify datasource and variables
    dataSource = setup.T;
    
    % If time series is longer than 365 days do bias correction
    % Mean temperature
    if sum(~isnan(statData.WMO.TEMP))>365
        var2BC = {'TEMP'};
        [T,slopeT,ViT] = biasCorrection(T, R_T, timeT', statData, var2BC, dataSource); % T is now bias corrected!
    end
    
    % Minimum temperature
    if sum(~isnan(statData.WMO.MIN))>365
        var2BC = {'MIN'};
        [maxT,slopeTmax,ViTmax] = biasCorrection(maxT, R_T, timeT', statData, var2BC, dataSource); % T is now bias corrected!
    end
    
    % Maximum temperature
     if sum(~isnan(statData.WMO.MAX))>365
        var2BC = {'MAX'};
        [minT,slopeTmin,ViTmin] = biasCorrection(minT, R_T, timeT', statData, var2BC, dataSource); % T is now bias corrected!
     end
        
    %Extract maximum 20 days
    if size(T,3)>20
        timeT=timeT(end-19:end);
        maxT=maxT(:,:,end-19:end);
        minT=minT(:,:,end-19:end);
        T=T(:,:,end-19:end);
    end
       
    
    
    
    %% 4. Compute potential ET from GDAS temperatures
    
    %calculate  centered latitude grid coordinates vector
    latT=zeros(size(T,1),size(T,2));
    for i=1:R_T.RasterSize(1,1)
        latT(i,:) = R_T.Latlim(1,2) - i*R_T.DeltaLat + R_T.DeltaLat/2;
    end
    
    %Trim Kc to GDAS
    [KcGDAS,R_Kc]=maptrims(Kc,R_NDVI,R_T.LatitudeLimits,R_T.LongitudeLimits);
    
    %calculate potential ET
    etG = single(potET(T,maxT,minT,latT,timeT')); %calculate pot. ET
    for l=1:size(etG,3)
        etGDAS(:,:,l)= imresize(etG(:,:,l),size(KcGDAS)).*KcGDAS;
    end
    R_ET=R_Kc;
       
   
    %save potential ET
    cd('../processed/temp')
    save(strcat(setup.T,'.mat'),'etGDAS', 'T','maxT','minT','R_ET','timeT', 'R_T')
    
    %clear variables
    clearvars  i maxT minT T latT etG etGDAS l R_ET R_T timeT timeNDVI KcGDAS R_Kc SM200 ET0
    
    
    %% 5. Bias correct GFS
    
    % load GFS data
    cd('../../raw')
    load(setup.F);
    
    % Specify datasource and variables
    dataSource = 'GFS';
    
    % If time series is longer than 365 days do bias correction
    % Mean temperature
    if sum(~isnan(statData.WMO.TEMP))>365
        var2BC = {'TEMP'}; 
        [F,slopeFT,ViFT] = biasCorrection(F, F{1}.R_F, size(F,2), statData, var2BC, dataSource); % T is now bias corrected
    end
    
    % Precipitation
    if sum(~isnan(statData.WMO.PRCP))>365
        var2BC = {'PRCP'}; 
        [F,slopeFP,ViFP] = biasCorrection(F, F{1}.R_F, size(F,2), statData, var2BC, dataSource); % T is now bias corrected!
    end
    
     % Mean temperature
    if sum(~isnan(statData.WMO.MIN))>365
        var2BC = {'MIN'}; 
        [F,slopeFminT,ViTmin] = biasCorrection(F, F{1}.R_F, size(F,2), statData, var2BC, dataSource); % T is now bias corrected!
    end
    
     % Mean temperature
    if sum(~isnan(statData.WMO.MAX))>365
        var2BC = {'MAX'}; 
        [F,slopeFmaxT,ViTmax] = biasCorrection(F, F{1}.R_F, size(F,2), statData, var2BC, dataSource); % T is now bias corrected!
    end
    
    %% 6.Compute potential ET from GFS temperatures
    
    if length(F)>5
    %extract last day
    F=F(end-4:end);
    end
    
    %centered latitude grid coordinates vector
    latF=zeros(size(F{1}.FT,1),size(F{1}.FT,2));
    for i=1:F{1}.R_F.RasterSize(1,1)
        latF(i,:) = F{1}.R_F.Latlim(1,2) - i*F{1}.R_F.DeltaLat + F{1}.R_F.DeltaLat/2;
    end
    
    %calculate pot. ET
    for l=1:length(F)
        d{l}.FET = single(potET(F{l}.FT,F{l}.FmaxT,F{l}.FminT,latF,F{l}.timeF));
    end
        
    
    %multiply with Kc
    for m=1:length(F)
        for l=1:size(F{m}.FT,3)
            for k=1:size(F{m}.FT,4)
                F{m}.FET(:,:,l,k)=imresize(d{m}.FET(:,:,l,k),size(Kc)).*single(Kc);
            end
        end
    end
    
    %adjust Georeference
     for l=1:length(F)
        F{l}.R_FP=F{l}.R_F;
        F{l}=rmfield(F{l},'R_F');
        F{l}.R_FET=R_NDVI;
     end
    
    %save potential ET
    cd('../processed/temp')
    save(strcat(setup.F,'.mat'),'F')
    
    %clear variables
    clearvars i latF FET F k l Kc R_NDVI d
    
    
    
    %% 7. Adjust, correct and downscale FEWS P to Subcatchments
    
    % read subcatchments
    SUB= setup.shp;
    
    % load FEWS
    cd('../../raw')
    load(strcat(setup.P,'.mat')); %load data
    
    %check if DEM is available
    cd('./dem')
    if exist(fullfile(cd, 'DEM_FEWS.mat'), 'file')==0
        % GET DEM FROM NASA WORLDWIND
        [DEM_FEWS,R_DEM_FEWS] = getDEM(R_P.LatitudeLimits, R_P.LongitudeLimits, 10);
        DEM_FEWS=double(DEM_FEWS);
        % Resize and sample
        [DEM_FEWS,R_DEM_FEWS] = resizem(DEM_FEWS,[size(P,1) size(P,2)],R_DEM_FEWS);
        % Save
        save DEM_FEWS.mat DEM_FEWS R_DEM_FEWS
    else
        load('DEM_FEWS.mat')
    end
    
    %check if nearest neighbor for FEWS raster has already been calculated
    if exist(fullfile(cd, 'iEiC_FEWS.mat'), 'file')==0
            [iE_FEWS, iC_FEWS] = findNeighbor(1500,R_P,DEM_FEWS);
            save iEiC_FEWS.mat iE_FEWS iC_FEWS;
    else
        load('iEiC_FEWS.mat');
    end
   
    % Adjust P to elevation
    % specify parameter
    param{1}.k1=0.4;
    param{1}.k2=0;
    [P] = adjMeteodata(P, param, 1, DEM_FEWS, iC_FEWS, iE_FEWS); % adjust to elevation  

    cd('..')
         
    
    % Bias correct
    dataSource = setup.P;
    
    % Precipitation %% Gives totally unrealistic values
    if sum(~isnan(statData.WMO.PRCP))>365
        var2BC = {'PRCP'};
        [P,slopeP,ViP] = biasCorrection(P, R_P, timeP', statData, var2BC, dataSource); % T is now bias corrected!
    end
   
    
    %Extract maximum 20 days
    if size(P,3)>20
        timeP=timeP(end-19:end);
        P=P(:,:,end-19:end);
    end
    newTFEWS=timeP;
    
     
    % downscale P to subcatchments
    cd('../processed/sub')
    if exist(fullfile(cd, 'SI_P.mat'), 'file')==0
        
       a_P = getSubIndex(SUB,P,R_P);
       save SI_P.mat a_P
        
    else
        load('SI_P.mat')
    end
    [new_P] = getSubValues(a_P,SUB,P,R_P);
    
    %Check if there is already a file
    if exist(fullfile(cd,'sub_P.mat'), 'file')==0
        Psub=new_P;
        timeP=newTFEWS;
        save sub_P.mat Psub timeP;
    else
        load('sub_P.mat')
        if newTFEWS(end)~= timeP(end)
           nNew=newTFEWS(end)-timeP(end); 
           Psub(end+1:end+nNew,:)=new_P(end-nNew+1:end,:);
           timeP(end+1:end+nNew)=newTFEWS(end-nNew+1:end);
        end
        save sub_P.mat Psub timeP;
    end
    
    clearvars sub_P timeP P R_P shape nNew  new_P a_P_FEWS newTFEWS
    
    %% Downscale GDAS to Subcatchments
    
    % load GDAS pET
    cd('../temp')
    load(strcat(setup.T,'.mat'))
    cd('../sub')
    
         
    %update date
    newTGDAS=timeT;
    
    % Potential Evapotranspiration
    if exist(fullfile(cd,'SI_pET.mat'), 'file')==0
        
        a_pET = getSubIndex(SUB,etGDAS,R_ET);
        save SI_pET.mat a_pET
        
    else
        load('SI_pET.mat')
    end
    [new_pET] = getSubValues(a_pET,SUB,etGDAS,R_ET);
    
    %Check if there is already a file
    if exist(fullfile(cd,'sub_pET.mat'), 'file')==0
        pETsub=new_pET;
        timeT=newTGDAS;
        save sub_pET.mat pETsub timeT;
    else
        load('sub_pET.mat')
        if newTGDAS(end)~= timeT(end)
           nNew=newTGDAS(end)-timeT(end); 
           pETsub(end+1:end+nNew,:)=new_pET(end-nNew+1:end,:);
           timeT(end+1:end+nNew)=newTGDAS(end-nNew+1:end);
        end
        save sub_pET.mat pETsub timeT;
    end
   
    
    % Temperatures
    if exist(fullfile(cd,'SI_T.mat'), 'file')==0
        
        a_T = getSubIndex(SUB,T,R_T);
        save SI_T.mat a_T
        
    else
        load('SI_T.mat')
    end
    
    [new_T] = getSubValues(a_T,SUB,T,R_T);
    [new_maxT] = getSubValues(a_T,SUB,maxT,R_T);
    [new_minT] = getSubValues(a_T,SUB,minT,R_T);
    
    %Check if there is already a file
    if exist(fullfile(cd,'sub_T.mat'), 'file')==0
        sub_T=new_T;
        sub_maxT=new_maxT;
        sub_minT=new_minT;
        timeT=newTGDAS;
        save sub_T.mat sub_T timeT;
        save sub_maxT.mat sub_maxT timeT;
        save sub_minT.mat sub_minT timeT;
    else
        load('sub_T.mat')
        load('sub_maxT.mat')
        load('sub_minT.mat')
        
        if newTGDAS(end)~= timeT(end)
           nNew=newTGDAS(end)-timeT(end); 
           sub_T(end+1:end+nNew,:)=new_T(end-nNew+1:end,:);
           sub_maxT(end+1:end+nNew,:)=new_maxT(end-nNew+1:end,:);
           sub_minT(end+1:end+nNew,:)=new_minT(end-nNew+1:end,:);
           timeT(end+1:end+nNew)=newTGDAS(end-nNew+1:end);
        end
        save sub_T.mat sub_T timeT;
        save sub_maxT.mat sub_maxT timeT;
        save sub_minT.mat sub_minT timeT;
    end
    
   
    clearvars sub_T sub_maxT sub_minT timeT R_T R_ET T a_T_GDAS a_pET_GDAS etGDAS m maxT minT newTGDAS new_T new_maxT new_minT new_pET sub_pET
    
   %% Downscale GFS to Subcatchments 
  
  
   % load GFS
    cd('../temp')
    load(strcat(setup.F,'.mat'))
    cd('../sub')
    
         
      
    % Temperatures and precipitation
    % Downscale new values
    if exist(fullfile(cd,'SI_P_F.mat'), 'file')==0
        
        a_P_F = getSubIndex(SUB,squeeze(F{1}.FP(:,:,1,1)),F{1}.R_FP);
        save SI_P_F.mat a_P_F
        
    else
        load('SI_P_F.mat')
    end
    for l=1:length(F)
        new{l}.FP = getSubValues(a_P_F,SUB,F{l}.FP,F{l}.R_FP); 
        new{l}.FT = getSubValues(a_P_F,SUB,F{l}.FT,F{l}.R_FP);
        new{l}.FmaxT = getSubValues(a_P_F,SUB,F{l}.FmaxT,F{l}.R_FP); 
        new{l}.FminT = getSubValues(a_P_F,SUB,F{l}.FminT,F{l}.R_FP);
        new{l}.timeF = F{l}.timeF;
    end
    % Update File
    if exist(fullfile(cd,'sub_FT.mat'), 'file')==0
        for l=1:length(new)
            FTsub{l}.FT=new{l}.FT;
            FTsub{l}.FmaxT=new{l}.FmaxT;
            FTsub{l}.FminT=new{l}.FminT;
            FTsub{l}.timeF=new{l}.timeF;
            FPsub{l}.FP=new{l}.FP;
            FPsub{l}.timeF=new{l}.timeF;
        end
        save sub_FP.mat FPsub;
        save sub_FT.mat FTsub;
       
    else
        load('sub_FP.mat')
        nNewFP=new{end}.timeF(1)-FPsub{end}.timeF(1);
        c=length(FPsub)+1;
            for m=length(new)-nNewFP+1:length(new)
                FPsub{c}.FP=new{m}.FP;
                FPsub{c}.timeF=new{m}.timeF;
                c=c+1;
            end
        load('sub_FT.mat')    
        nNewFT=new{end}.timeF(1)-FTsub{end}.timeF(1);
        c=length(FTsub)+1;
            for m=length(new)-nNewFT+1:length(new)
                FTsub{c}.FT=new{m}.FT;
                FTsub{c}.FmaxT=new{m}.FmaxT;
                FTsub{c}.FminT=new{m}.FminT;
                FTsub{c}.timeF=new{m}.timeF;
                c=c+1;
            end
        save sub_FP.mat FPsub;
        save sub_FT.mat FTsub;
       
    end
    
       
        
    % Potential Evapotranspiration
    % Downscale new values
    if exist(fullfile(cd,'SI_pET_F.mat'), 'file')==0
        
        a_pET_F = getSubIndex(SUB,squeeze(F{1}.FET(:,:,1,1)),F{1}.R_FET);
        save SI_pET_F.mat a_pET_F
        
    else
        load('SI_pET_F.mat')
    end
    for l=1:length(F)
        new{l}.FpET = getSubValues(a_pET_F,SUB,F{l}.FET,F{l}.R_FET);
    end
    % Update File
    if exist(fullfile(cd,'sub_FpET.mat'), 'file')==0
        for l=1:length(new)
            FpETsub{l}.FpET=new{l}.FpET;
            FpETsub{l}.timeF=new{l}.timeF;
        end
        save sub_FpET.mat FpETsub;
    else
        load('sub_FpET.mat')
        nNewFpET=new{end}.timeF(1)-FpETsub{end}.timeF(1);
        c=length(FpETsub)+1;
            for m=length(new)-nNewFpET+1:length(new)
                FpETsub{c}.FpET=new{m}.FpET;
                FpETsub{c}.timeF=new{m}.timeF;
                c=c+1;
            end
         save sub_FpET.mat FpETsub;   
    end
    
    
    %%
%     disp(datestr(today))
    clearvars today
    toc
 catch
     disp('Error occured, next try in 3 hours')
end











