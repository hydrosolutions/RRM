%% SOIL WATER STRESS INDEX (SWASI) - and thence, there was a name! :)

% Operational tool for calculating farmer (site)-specific soil water
% conditions. Requires carefully server-side implementation for 'online'
% calculations and updating (assimilation).
%
% SMS Trigger: STRESSINDEX siteID crop type soil type date (default=today)
%
% IMPORTANT QUESTIONS BELOW:
% > Thinking required with regard to triggering of tool. I.e. does the farmer
% automatically start an instance of SWASI or does he use a specific
% ICT channel trigger.
% > Assuming now that farmer-specific instances of SWASI are running on the
% server, then there should also be a way to stop/ reset a session.
% > are current states at t stored in the database and then querried for
% the t+1 calculation or the assimilation. this would require storage of
% the full ensemble!
% > How to deal with situations where a farmer has multiple fields. Would
% he then be running several instances? 
% > does the farmer only get back the SWASI if he sends new data or is
% there an independent way to query current conditions, i.e. remotely.
% > Reporting back one our of 5 classes (current understanding). Would
% maybe a class division as in the case of the soft soil moisture sensing
% make more sensing (not remembering why we decided to go for reporting
% back classes and not, say, ensemble mean states.
% > WHAT ARE MEANINGFUL INITAL PARAMETER STATES? 
% > CAN I (AS FARMER) BE 'STARTING' MY
% SWASI ANY TIME OR DOES THE GUY DEPEND ON RUNNING FROM THE VERY BEGINNING
% OF THE AG SEASON? PUT DIFFERENTLY, HOW DO I GET THE INITIAL CONDITIONS.
% GIVEN THE CURRENT SMS CHANNEL COMMUNICATION POSSIBILITIES, HOW CAN A) FIELD
% CAPACITY, B) WILTING POINT, C) ROOT DEPTH BE ENTERED (ASSUMING THAT DAY
% OF YEAR CAN BE CALCULATED AND THAT THE DEPLETION FACTOR FOR MAIZE IS
% STORED SOMEWHERE IN THE DATABASE (IS THAT SO? IF YES, WHERE?).
% > OTHER QUESTIONS RELATED TO HOW MET STATION DATA CAN BE ACCOUNTED FOR
% (BUT LETS DO THAT IN A VERSION 2).
%
% Comments: tobias siegfried, 30/04/2014
% Programming: Sebastian Stoll

%% Calculate Water Soil Balance stepwise

clear all, clc, close all

%% Load Climate Data @SEBI: WHICH SAMPLE YEAR IS THIS?
% @SEBI: WOULD THESE BE PICKED FROM THE DATABASE?

if ismac == 1
    cd('~/Dropbox (hydrosolutions)/RS_CA_iMoMoBridgingPhase/waterStress/')
else
    cd('C:\Users\Sebastian\Dropbox (hydrosolutions)\iMoMoBridgingPhase\waterStress')
end

load('testData')
% @SEBI: What kind of data is this and where does it come from? Do these
% data need to be updated regularly.

%% Define Parameters

noP = 20; % Number of parameter Ensembles
noC = 20; % Number of climate Ensembles
actID = 80; % Day of the year
p = 0.55; % Depletion factor for Maize - @SEBI: WHAT ABOUT FOR OTHER CROPS?
qFC = 0.4; % WC at field capacity
qWP = 0.1; % WC at wilting point
R = zeros(noC,noP); % Runoff
D = 0.7; % Root depth
startWC = 0.5; % fraction of WC filling at beginning
% TAW=1000*(qFC-qWP)*D; % total available water in mm
% RAW=p*TAW; %RAW=Readily available water,  p=depletion factor
% Dr=RAW*startWC; %Root zone depletion
% Irr=zeros(100); %Irrigation Demand

%% Generate parameter Ensembles (@SEBI: PLEASE DESCRIBE EACH AND EVERY VARIABLE!)

Ep1 = normrnd(p,0.05,1,noP);
EqFC1 = normrnd(qFC,0.02,1,noP);
EqWP1 = normrnd(qWP,0.02,1,noP);
EqWP1(EqWP1 <= 0.02) = 0.02;
ED1 = normrnd(D,0.05,1,noP);
EstartWC1 = normrnd(startWC,0.05,1,noP);
TAW1 = 1000 .* (EqFC1 - EqWP1) .* ED1; % total available water in mm
MaxWC1 = 1000 .* EqFC1 .* ED1;
RAW1 = Ep1 .* TAW1; %RAW=Readily available water,  p=depletion factor
Dr1 = RAW1 .* EstartWC1; %Root zone depletion
Irr = zeros(noC,noP); %Irrigation Demand

for l=1:noP
    Ep(l,:) = Ep1;
    EqFC(l,:) = EqFC1;
    EqWP(l,:) = EqWP1;
    ED(l,:) = ED1;
    EstartWC(l,:) = EstartWC1;
    TAW(l,:) = TAW1;
    MaxWC(l,:) = MaxWC1;
    RAW(l,:) = RAW1;
    Dr(l,:) = Dr1;
end

IR(:,:,1) = Irr;
DR(:,:,1) = Dr;
WC(:,:,1) = (MaxWC-Dr) ./ MaxWC .* EqFC;

%% Generate Climate Ensembles @SEBI: WHAT IS REQUIRED IN TERMS OF KEEPING UPDATES OF THESE STATES?
% JUST QUERYING THE DATABASE FOR REMOTELY SENSED INFORMATION? WHAT ABOUT
% DATA FROM NEARBY MET STATIONS?

ETc = ETC{actID}(1,:);
Pc = PC{actID}(1,:);

for l=1:noC
    EETc1(l,:)=normrnd(ETc,0.5);
end

for l=1:noC
    EPc1(l,:)=normrnd(Pc,2.5);
end

EPc1(EPc1<0)=0;

for k=1:size(ETc,2)
    
    for l=1:noC
        
        EETc(:,l,k)=EETc1(:,k);
        EPc(:,l,k)=EPc1(:,k);
        
    end
    
end

%% Run Model (@SEBI: PLEASE DESCRIBE VARIABLES)

for l=1:size(EPc,3)
    [Dr,Irr]=StepSoilB(TAW,Ep,RAW,Dr,R,EPc(:,:,l),EETc(:,:,l),Irr);
    DR(:,:,l+1)=Dr;
    IR(:,:,l+1)=Irr;
    WC(:,:,l+1)=(MaxWC-Dr)./MaxWC.*EqFC;
end

% % Plot water contents
% for k=1:size(WC,1)
%     hold on
%     for l=1:size(WC,2)
%         plot(squeeze(WC(k,l,:)))
%     end
% end

% Make 'fake' observations
x_true1 = squeeze(mean(mean(WC,1)));

%% Generate Parameter Ensembles (@SEBI VARIABLE)

Ep1 = normrnd(p,0.05,1,noP);
EqFC1 = normrnd(qFC,0.02,1,noP);
EqWP1 = normrnd(qWP,0.02,1,noP);
EqWP1(EqWP1 <= 0.02) = 0.02;
ED1 = normrnd(D,0.05,1,noP);
EstartWC1 = normrnd(startWC,0.05,1,noP);
TAW1 = 1000 .* (EqFC1 - EqWP1) .* ED1; % total available water in mm
MaxWC1 = 1000 .* EqFC1 .* ED1;
RAW1 = Ep1 .* TAW1; %RAW=Readily available water,  p=depletion factor
Dr1 = RAW1 .* EstartWC1; %Root zone depletion
Irr = zeros(noC,noP); %Irrigation Demand

for l = 1 : noP
    
    Ep(l,:) = Ep1;
    EqFC(l,:) = EqFC1;
    EqWP(l,:) = EqWP1;
    ED(l,:) = ED1;
    EstartWC(l,:) = EstartWC1;
    TAW(l,:) = TAW1;
    MaxWC(l,:) = MaxWC1;
    RAW(l,:) = RAW1;
    Dr(l,:) = Dr1;
    
end

IR(:,:,1) = Irr;
DR(:,:,1) = Dr;

%% Data Assimilation (@SEBI: PLEASE DOCUMENT, INCLUDING ALL VARIABLES)

m = noP * noC; %number of ensembles
r = 0.005;
infl = 1.35;
upstep = 14;
id(1) = 0;
cc =1 ; % counter

for step = 1 : 125 % time-steps (daily) - @SEBI: Why 125?
    
    % Propagate the ensemble, calculate ensemble observations.
    HE = [];
    [Dr,Irr] = StepSoilB(TAW,Ep,RAW,Dr,R,EPc(:,:,step),EETc(:,:,step),Irr);
    Wc = (MaxWC - Dr) ./ MaxWC .* EqFC;
    DR(:,:,step+1) = Dr;
    IR(:,:,step+1) = Irr;
    WC(:,:,step+1) = (MaxWC-Dr) ./ MaxWC .* EqFC;
    id(step+1) = 0;
    
    % ASSIMILATION STEP
    if mod(step, upstep) == 0 %if assimilation step
        
        HE=Wc(:)'; %ensemble states at current timestep
        y = SoilB_observe(x_true1(step), r); %produce observation from true field and stdv
        
        % calculate ensemble observation anomalies
        Hx = mean(HE,2); % Ensemble Mean
        HA = HE - Hx; % Deviation of Ensembles from Ensemble Mean
        dy = y - Hx; %Deviation of Ensemble mean from observation
        
        % calculate standardised innovation and standardised ensemble anomalies
        s = dy / (r * sqrt(m - 1)); %Deviation of Ensemble mean from observation
        S = HA / (r * sqrt(m - 1)); %Deviation of Ensembles from Ensemble Mean
        [U, L] = svd(speye(m) + S' * S);
        l = diag(L);
        
        % note: U and l will be reused for calculating T and inv(T)
        G = U * diag(1 ./ l) * U';
        b = G * S' * s;
        dx = HA * b;
        
        x2 = Hx; % Ensemble mean
        A2 = HE - repmat(Hx, 1, m);% Deviation from Ensemble Mean
        x2 = x2 + A2 * b;
        T = sqrtm(G);
        A2 = A2 * (T * infl);
        E2 = repmat(x2, 1, m) + A2; %update to new states
        
        %convert to root zone depletion
        Wc = (MaxWC-Dr) ./ MaxWC .* EqFC;
        WcE = reshape(E2,noC,noP);
        Dr = (WcE ./ EqFC .* MaxWC - MaxWC) * -1;
        
        %IDs
        id(step+1) = 1;
        stepId(cc) = step + 1;
        cc = cc + 1;
        
    end
    
end

id=logical(id);

% Plot water contents
for k = 1 : size(WC,1)
    
    figure(2)
    hold on
    for l = 1 : size(WC,2)
        h1 = plot(squeeze(WC(k,l,:)));
    end
    grid minor
    
end

h2 = plot(stepId+1,x_true1(id'),'rx','MarkerSize',15,'linewidth',2);
legend([h1 h2],'model ensemble','observation')
xlabel('days since planting');
ylabel('soil moisture (VWC)');

%% USED BY SEBI FOR HIS PRESENTATION FIGURES
% for o=1:7:length(x_true1)
%     f1=figure(4);
%
%     subplot(3,1,1)
%         h1=plot(ETc(1:o),'r','linewidth',2);
% %         xlabel('days since planting');
%         ylabel('crop ET (mm)');
%         ylim([0 7])
%         xlim([0 130])
%         title(strcat('Day',{' '},num2str(o)))
%         grid minor
%
%     subplot(3,1,2)
%         h2=b(Pc(1:o),'b');
% %         xlabel('days since planting');
%         ylabel('P (mm)');
%         ylim([0 60])
%         xlim([0 130])
% %         title(strcat('Day',{' '},num2str(o)))
%          grid minor
%
%     subplot(3,1,3)
%         h3=plot(x_true1(1:o),'g','linewidth',2);
%         xlabel('days since planting');
%         ylabel('soil moisture (VWC)');
%         ylim([0.1 0.5])
%         xlim([0 130])
% %         title(strcat('Day',{' '},num2str(o)))
%         grid minor
%
%
%     saveas(f1, strcat('Single',num2str(o)),'png');
% end

































































































































