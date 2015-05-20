function [prm] = model_setprm(prm)
%% model_setprm sets additional parameters for the basin specific custmprm
% set
% tobias siegfried, 27//11/2014

prm.customprm.F = (prm.customprm.Fmin + prm.customprm.Fmax) / 2;
prm.nprm = prm.n - prm.p; % assuming that these are the number of parameters
prm.periodic = 1;

% rrModel-iMoMo add-on!
load(['../resources/samples/' prm.customprm.fname_samples],'time');
load(['../resources/samples/' prm.customprm.fname_samples],'P');
load(['../resources/samples/' prm.customprm.fname_samples],'ET'); % ET actual
load(['../resources/samples/' prm.customprm.fname_samples],'Tmin');
load(['../resources/samples/' prm.customprm.fname_samples],'Tmax');
load(['../resources/samples/' prm.customprm.fname_samples],'Tmean');
load(['../resources/samples/' prm.customprm.fname_samples],'Lat');
load(['../resources/samples/' prm.customprm.fname_samples],'SM'); % soil moisture
load(['../resources/samples/' prm.customprm.fname_samples],'cID'); % id of active catchments over which assim. takes place


prm.customprm.time = time;
prm.customprm.P = P;
prm.customprm.ET = ET;
prm.customprm.Tmin = Tmin;
prm.customprm.Tmax = Tmax;
prm.customprm.cID = cID;

if prm.customprm.snow
    load(['./' prm.customprm.fname_samples],'SN'); % snow
    prm.customprm.SN = SN;
end
    
try
%     temp = load(prm.customprm.geometry);
    temp = shaperead(prm.customprm.geometry);
    prm.shp = temp;
catch
    disp('ERROR loading geometry')
end

try
    prm.customprm.Tmean = Tmean;
catch
    prm.customprm.Tmean = 0;
end

try
    prm.customprm.SM = SM;
catch
    prm.customprm.Tmean = 0;
end

try
    prm.customprm.Q = Q;
catch
    prm.customprm.Q = zeros(size(P));
end

try
    prm.customprm.G = G;
catch
    prm.customprm.G = zeros(size(P));
end

prm.customprm.Lat = Lat;


try
    load([prm.customprm.loc_connMat],'connMat');
    load([prm.customprm.loc_connMat],'connMatInd');
    prm.customprm.connMat = connMat;
    prm.customprm.connMatInd = connMatInd;
catch
    disp('Connection Matrix not yet available!')
end

return
end