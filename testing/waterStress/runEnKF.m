%% Calculate Water Soil Balance stepwise

clear all
%% load climate data
if ismac == 1
cd('C:\Users\Jules\Dropbox\Dokumente\Jobs\hydrosolutions\RRM\testing\waterStress')
else
cd('C:\Users\Jules\Dropbox\Dokumente\Jobs\hydrosolutions\RRM\testing\waterStress')
end

load('testData')

%% Define Parameter
noP=20; %number of parameter Ensembles
noC=20; %number of climate Ensembles


%Soil water balance model parameters
actID=80; %day of the year 
p=0.55; %Depletion factor  for Maize
qFC=0.4; %WC at field capacity 
qWP=0.1; %WC at wilting point
R=zeros(noC,noP); %Runoff
D=0.7; %Root depth
startWC=1; %fraction of WC filling at beginning
% TAW=1000*(qFC-qWP)*D; % total available water in mm
% RAW=p*TAW; %RAW=Readily available water,  p=depletion factor
% Dr=RAW*startWC; %Root zone depletion
% Irr=zeros(100); %Irrigation Demand

%% Generate parameter Ensembles
Ep1 = normrnd(p,0.05,1,noP);
EqFC1 = normrnd(qFC,0.02,1,noP);
EqWP1 = normrnd(qWP,0.02,1,noP);
    EqWP1(EqWP1<=0.02)=0.02;  
ED1 = normrnd(D,0.05,1,noP);
EstartWC1=normrnd(startWC,0.05,1,noP);
TAW1=1000.*(EqFC1-EqWP1).*ED1; % total available water in mm
MaxWC1=1000.*EqFC1.*ED1;
RAW1=Ep1.*TAW1; %RAW=Readily available water,  p=depletion factor
Dr1=RAW1.*EstartWC1; %Root zone depletion
Irr=zeros(noC,noP); %Irrigation Demand
for l=1:noP
    Ep(l,:)=Ep1;
    EqFC(l,:)=EqFC1;
    EqWP(l,:)=EqWP1;
    ED(l,:)=ED1;
    EstartWC(l,:)=EstartWC1;
    TAW(l,:)=TAW1;
    MaxWC(l,:)=MaxWC1;
    RAW(l,:)=RAW1;  
    Dr(l,:)=Dr1;
end
IR(:,:,1)=Irr;
DR(:,:,1)=Dr;
WC(:,:,1)=(MaxWC-Dr)./MaxWC.*EqFC;
%% Generate climate Ensembles
ETc=ETC{actID}(1,:);
Pc=PC{actID}(1,:);

for l=1:noC
EETc1(l,:)=normrnd(ETc,1.5);
end

for l=1:noC
EPc1(l,:)=normrnd(Pc,3.5);
end
EPc1(EPc1<0)=0;

for k=1:size(ETc,2)
    for l=1:noC
        EETc(:,l,k)=EETc1(:,k);
        EPc(:,l,k)=EPc1(:,k);
    end
end


%% Run model to generate X_true
for l=1:size(EPc,3)
    [Dr,Irr]=StepSoilB(TAW,Ep,RAW,Dr,R,EPc(:,:,l),EETc(:,:,l),Irr);
    DR(:,:,l+1)=Dr;
    IR(:,:,l+1)=Irr;
    WC(:,:,l+1)=(MaxWC-Dr)./MaxWC.*EqFC;
    
end


% % Plot water contents
for k=1:size(WC,1)
    figure(2);
    hold on
    for l=1:size(WC,2)
        h1=plot(squeeze(mean(mean(WC,1))));
    end
end

%make fake observations
x_true1=squeeze(mean(mean(WC,1)));
%% Generate parameter Ensembles
Ep1 = normrnd(p,0.05,1,noP);
EqFC1 = normrnd(qFC,0.02,1,noP);
EqWP1 = normrnd(qWP,0.02,1,noP);
    EqWP1(EqWP1<=0.02)=0.02;  
ED1 = normrnd(D,0.05,1,noP);
EstartWC1=normrnd(startWC,0.05,1,noP);
TAW1=1000.*(EqFC1-EqWP1).*ED1; % total available water in mm
MaxWC1=1000.*EqFC1.*ED1;
RAW1=Ep1.*TAW1; %RAW=Readily available water,  p=depletion factor
Dr1=RAW1.*EstartWC1; %Root zone depletion
Irr=zeros(noC,noP); %Irrigation Demand
for l=1:noP
    Ep(l,:)=Ep1;
    EqFC(l,:)=EqFC1;
    EqWP(l,:)=EqWP1;
    ED(l,:)=ED1;
    EstartWC(l,:)=EstartWC1;
    TAW(l,:)=TAW1;
    MaxWC(l,:)=MaxWC1;
    RAW(l,:)=RAW1;  
    Dr(l,:)=Dr1;
end
IR(:,:,1)=Irr;
DR(:,:,1)=Dr;



%% Data assimilation

m=noP*noC;
r=0.01;
infl=1.35;
upstep=21;
id(1)=0;
cc=1; %counter

for step = 1 : 125
    
        
    % propagate the ensemble, calculate ensemble observations
    HE = [];
    [Dr,Irr]=StepSoilB(TAW,Ep,RAW,Dr,R,EPc(:,:,step),EETc(:,:,step),Irr);
    Wc=(MaxWC-Dr)./MaxWC.*EqFC;
    DR(:,:,step+1)=Dr;
    IR(:,:,step+1)=Irr;
    WC(:,:,step+1)=(MaxWC-Dr)./MaxWC.*EqFC;
    id(step+1)=0;
    WCbefore=WC;
    
    
   
    
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
            Wc=(MaxWC-Dr)./MaxWC.*EqFC;
            WcE= reshape(E2,noC,noP);
            Dr=(WcE./EqFC.*MaxWC-MaxWC)*-1;
            
            %IDs
            id(step+1)=1;
            stepId(cc)=step+1;
            cc=cc+1;
        end
       
end

id=logical(id);
c=1;
% Plot water contents
for k=1:size(WC,1)
           
    figure(2)
    hold on
    for l=1:size(WC,2)
        h4=plot(squeeze(mean(mean(WCbefore,1))));
        all(:,c)=squeeze(WC(k,l,:));
        c=c+1;
   end
    grid minor
end
h2=plot(mean(all,2),'g','linewidth',2);
h3=errorbar(stepId+1,x_true1(id'),ones(1,length(x_true1(id')))*r,'rx','MarkerSize',15,'linewidth',2);
legend([h1 h2 h3 h4],'true state','ensemble mean','observation','before')
xlabel('days since planting');
ylabel('soil moisture (VWC)');


%%
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


















    
    
    












































































































