function [x] = budykoStepFSVec(prm,x,p,true_field)
% Usage:
%                   [x,varargout] = budykoStepFS(prm,x,p,true_field)
%
% Input:            x: state variables at time t, {q,S,G} ? NOW: size = p x m
%                   p: parameter values at time t ? NOW: size = p x m
%
% Output:           x(t+1)
%
% Created:          15 / 09 / 2011
%
% Last modified:    12 / 11 / 2012
%
% Author:           Tobias Siegfried
%
% Purpose:          Performs a one timestep water balance calculation
%                   using the Budyko scheme as reported in Zhang, 2008.
%
% Revisions:        29 / 04 / 2013, removed the prm.themiSUB loading and
%                   added to get_prm file.
%                   01 / 02 / 2013, subcatchment flow routing now included.
%                   13 / 12 / 2012, vectorize function for large speed gains.
%                   12 / 12 / 2012 switch (0/1) for snow and restoring the states for coherence.
%                   (work started).
%                   11 / 12 / 2012: Correct the date entries for grabbing
%                   the model drivers (P and T).
%                   04 / 12 / 2012: Tmean is available (I/O) - use it if available!
%                   13 / 11 / 2012: debugging ok.
%                   06 / 11 / 2012: adaptations for iMoMo, including snow
%                   storage plus the recognition of the fact that ETact is also a state (in
%                   light of the fact that we will be assimilating remote sensing data into
%                   the model. So, from a 3 state model, we will go to a 5
%                   state model (x = [Q;S(t);G(t);SR(t);ET], plus there are two additional parameters
%                   (now 6 in total). p(4) is rain snow transition temperature an par(5) is
%                   the melting factor. Function is now called budykoStepFS, where FS stands for full states.
%                   06 / 11 / 2011: addressed the initial condition problem, numerous other fixes
%                   18 / 10 / 2011: returns the ET(t) as an optional
%                   argument - obsolete.
%% date (day of the year) - note format: 'yyyymmdd'
dn = prm.current.doy;
nC = prm.customprm.nC;
cID = prm.customprm.cID;
sizX2 = size(x,2);
%% atmospheric drivers
tFromBeg = prm.current.serialDate - prm.customprm.serialDateStart + 1;
P = prm.customprm.P(tFromBeg,:)'; % precipitation, at time t
P = repmat(P,1,sizX2);


if prm.customprm.calculateET
    Tmin = prm.customprm.Tmin(tFromBeg,:)'; % Temperature, at time t
    Tmax = prm.customprm.Tmax(tFromBeg,:)'; % Temperature, at time t
    if prm.customprm.TmeanAvail
        Tmean = prm.customprm.Tmean(tFromBeg,:)';
    else
        Tmean = (Tmin+Tmax) / 2; % mean temperature at time t
    end
end

t = 2; % only operates on the {t,t-1} steps where t=2 (=prm.current.step)
S = zeros(2 * nC, sizX2); % unsaturated zone storage
G = zeros(2 * nC, sizX2); % groundwater storage
SN = zeros(2 * nC, sizX2); % snow storage
%% initial condition?
if prm.current.step == 1 % this obv. needs to be improved under the restarting condition.
    if strcmp(prm.customprm.runType,'training')
        S = repmat(prm.customprm.S0(:), 2, sizX2); % [t; t-1]
        G = repmat(prm.customprm.G0(:), 2, sizX2); % [t; t-1]
        if prm.customprm.snow
            SN = repmat(prm.customprm.SN0, 2 * nC, sizX2); % [t; t-1]
        end
    elseif strcmp(prm.customprm.runType,'online')
        S(1:2:end-1,:) = prm.storage.S0G0(1:2:end-1,:);
        G(1:2:end-1,:) = prm.storage.S0G0(2:2:end,:);
    else
        error('Wrong run type. Please choose either ''training'' or ''online''!');
    end
    
else
    % Initial conditions for the storage compartments
    if ~true_field
        if ~prm.customprm.snow
            S(t-1:2:end-1,:) = squeeze(prm.storage.SG(1:2:end-1,prm.current.step-1,:)); %S
            G(t-1:2:end-1,:) = squeeze(prm.storage.SG(2:2:end,prm.current.step-1,:));   %G
        elseif prm.customprm.snow
            S(t-1:3:end-1,:) = squeeze(prm.storage.SG(1:3:end-2,prm.current.step-1,:)); %S
            G(t-1:3:end-1,:) = squeeze(prm.storage.SG(2:3:end-1,prm.current.step-1,:)); %G
            SN(t-1:3:end-1,:) = prm.storage.SG(3:3:end,prm.current.step-1,:);           %SN
        end
    else
        S(t-1:2:end-1,:) = prm.trueX(2:8:end,tFromBeg-1); % check
        G(t-1:2:end-1,:) = prm.trueX(3:8:end,tFromBeg-1);
        if prm.customprm.snow
            SN(t-1:2:end-1,:) = prm.trueX(4:prm.n:end,tFromBeg-1);
        end
    end
end
%% Potential ET (= E0) at the particular day/date - no vec needed!
if prm.customprm.calculateET
    E0 = hargreavesPET(Tmax,Tmin,Tmean,prm.customprm.Lat,dn);
else
    E0 = prm.customprm.ET';
end
E0 = repmat(E0,1,sizX2); % for multiple catchments
%% Temperature & Snow Melt (ATTN.: NOT YET VECTORIZED and not yet fit for multiple catchment assimilation!)
if prm.customprm.snow
    if Tmean <= p(5)    % negative temperature  -> all precipitation for subbasin stored as snow
        SN(t:2:end,:) = SN(t-1:2:end-1,:) + P; P = 0;
    else                % snow of subbasin melts
        m = p(6) * (Tmean - p(5)); % potential snowmelt
        m = min(m,SN(t-1)); % max snow melt = snow storage from t-1
        P = P + m;
        SN(t) = SN(t-1) - m;  if m < 0, disp('ERROR: m < 0'), return, end
    end
end
%% Catchment Rainfall Retention X
cRet = P <= prm.customprm.limCutoff;
X = zeros(nC, sizX2);

p1 = p(1:(prm.n-prm.nprm):end-3,:);
p4 = p(4:(prm.n-prm.nprm):end,:);
ST1 = S(t-1:2:end-1,:);
if any(~cRet)
    A=p4(~cRet) - ST1(~cRet) + E0(~cRet);
    B= P(~cRet);
    alph=p1(~cRet);
        
    w = 1 ./ (1 - alph);
    C = (A)./(B);
    te = C .^ w;
    if any(te(:) > 10^10) || any(isnan(te)) % XXX: stop hard-wiring and adress this properly in the parameter definition files
       f = C; f(C>=1) = 1;
       f=double(f);
    else
       f = 1 + C - (1 + te) .^ (1./w); f(f>=1)=1;
       f=double(f);
    end
    X(~cRet) = P(~cRet) .* f;
end

%% Direct Runoff
Qd = P - X;
%% Water Availability
W = X + S(t-1:2:end-1,:);
Y = zeros(nC,sizX2); ET = Y; % matrix prep.

watAv = W <= prm.customprm.limCutoff;
% Y and ET set = 0 already. just deal with the ~0 case
p2 = p(2:(prm.n-prm.nprm):end-2,:);
if any(~watAv)
      
    A=E0(~watAv)+ p4(~watAv);
    B= W(~watAv);
    alph=p2(~watAv);
        
    w = 1 ./ (1 - alph);
    C = (A)./(B);
    te = C .^ w;
    if any(te(:) > 10^10) || any(isnan(te)) % XXX: stop hard-wiring and adress this properly in the parameter definition files
       f = C; f(C>=1) = 1;
       f=double(f);
    else
       f = 1 + C - (1 + te) .^ (1./w); f(f>=1)=1;
       f=double(f);
    end
    Y(~watAv) = W(~watAv) .* f;
    
    
    A=E0(~watAv);
    B= W(~watAv);
    alph=p2(~watAv);
        
    w = 1 ./ (1 - alph);
    C = (A)./(B);
    te = C .^ w;
    if any(te(:) > 10^10) || any(isnan(te)) % XXX: stop hard-wiring and adress this properly in the parameter definition files
       f = C; f(C>=1) = 1;
       f=double(f);
    else
       f = 1 + C - (1 + te) .^ (1./w); f(f>=1)=1;
       f=double(f);
    end
    ET(~watAv) = W(~watAv)  .* f;
end

%% Rootzone Storage
S(t:2:end,:) = Y - ET;
S(S<0)=0; %Catch numeric errors
%% Groundwater Recharge
R = W - Y;
%% Groundwater Balance
G(t:2:end,:) = (1 - p(3:(prm.n-prm.p):end-1,:)) .* G(t-1:2:end-1,:) + R;
%% Baseflow (lin. reservoir)
Qb = p(3:(prm.n-prm.p):end-1,:) .* G(t-1:2:end-1,:);
%% Total Combined Runoff
Q = Qb + Qd;

% convert to m3/s
for c = 1 : prm.customprm.nC
    qConvFac(c) = 1/(1000*24*3600) * prm.shp(prm.customprm.cID(c)).Area;
    Qq(c,:)=Q(c,:)*qConvFac(c);
end

%% Total Routed Combined Runoff
QRout = Qq;
for i = 1 : size(prm.customprm.connMatInd,1) % note: vectorization not possible - unfortunately!
    QRout(prm.customprm.connMatInd(i,2),:) = QRout(prm.customprm.connMatInd(i,2),:) + QRout(prm.customprm.connMatInd(i,1),:);
end
Q = QRout;
%% Final Output of State Variables
if prm.customprm.snow % :> not yet fit for multiple catchments...
    x = [Q;S(t,:);G(t,:);SN(t,:);ET]; % this now the new full state vector with snow.
else % no snow
    tQ = Q(:)'; tS = S(t:2:end,:); tS = tS(:)';
    tG = G(t:2:end,:); tG = tG(:)'; tET = ET(:)';
    x = [tQ;tS;tG;tET];
    x = reshape(x, nC * (prm.n-prm.p), sizX2);
end
%% Physics IO?
if ~isreal(x(:)) || any(isinf(x(:))) || any(isnan(x(:))) || any(x(:) < 0)
    [i,j] = find(x < 0);
    [i j]
    disp('Problem with Model Physics!? CHECK!'); keyboard
end
%% NESTED FUNCTIONS
%     function f = budykoCurveBIn(A,B,alph)
%         % See: Zhang et al. 2008 for explanation of the curve.
%         % REVISION: 11 / 10 / 2011: account for the odd cases and make it
%         % numerically more stable.
%         % tobias siegfried, 19 / 09 / 2011
%         
%         w = 1 ./ (1 - alph);
%         C = (A)./(B);
%         te = C .^ w;
%         if any(te(:) > 10^10) || any(isnan(te)) % XXX: stop hard-wiring and adress this properly in the parameter definition files
%             f = C; f(C>=1) = 1;
%         else
%             f = 1 + C - (1 + te) .^ (1./w);
%         end
%         
%     end
    function [Eo,lEo] = hargreavesPET(Tmax,Tmin,Tmean,Latitude,dn)
        % Potential evapotranspiration [mm/day], orig. from PBG code.
        T = [-20 -10 0 10 20 30 40];
        l = [2.549 2.525 2.501 2.477 2.453 2.430 2.406];
        %l_Tavg = interp1(T,l,Tmean);
        l_Tavg = qinterp1(T,l,Tmean);
        phi = Latitude'*2*pi/360;
        delta = SolarDeclination(dn); % solar declination in radians
        H0 = ETrad(delta,phi,dn);
        lEo = (0.0023.*H0.*(Tmax-Tmin).^0.5).*(Tmean+17.8);
        Eo = lEo ./ l_Tavg;
    end
    function delta = SolarDeclination(dn)
        %Solar declination [rad], orig from PBG
        delta = asin(0.4*sin((2*pi/365)*(dn-82)));
    end
    function H0 = ETrad(delta,phi,dn)
        %Extraterrestrial radiation [??]
        omega = 0.2618; %angular velocity of the earth's rotation [rad/hr]
        E0 = EccCorr(dn);
        T_SR = HourOfSunrise(delta,phi);
        H0 = 37.59*E0.*(omega*T_SR.*(sin(delta)*sin(phi))...
            +(cos(delta)*cos(phi).*sin(omega*T_SR)));
    end
    function T_SR = HourOfSunrise(delta,phi)
        %Hour of sunrise
        omega = 0.2618; %angular velocity of the earth's totation [rad/hr]
        T_SR = acos(-tan(delta).*tan(phi))./omega;
    end
    function E0 = EccCorr(dn)
        %Eccentricity correction of the earth's orbit
        E0 = 1+0.033*cos(2*pi*dn/365);
    end
end