function [Dr,Irr]=StepSoilB(TAW,p,RAW,Dr,R,Pc,ETc,Irr)
% Tobi: Sebi, please document your code in detail, otherwise not even the
% virgin mother marry, will ever find out what is what... this is
% especially crucial as Maxat/Diego should be working with this very soon!

Ks=(TAW-Dr)./((1-p).*TAW); % transpiration reduction factor dependent on available WC [0-1]
Ks(Dr<RAW)=1; % set to 1 if Dr<RAW

DP=Pc-R-(ETc.*Ks)-Dr; %Deep percolation
DP(Dr>=0)=0; % set to 0 if Dr>0

Dr=Dr-Pc-R+(ETc.*Ks)+DP;
DP(Dr<0)=Dr(Dr<0);
Dr(Dr<0)=0;
Irr=Dr-RAW;
Irr(RAW>Dr)=0;

return