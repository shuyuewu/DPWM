% Copyright (C) 2024 Shuyue Wu
% reference
% Shuyue Wu,Jianshi Zhao, Murugesu Sivapalan. A Parsimonious Daily Water Balance Model based on the Proportionality Hypothesis.(submitted to Journal of Hydrology)

load('rainplusmelt');
load('PET');
load('parameters');

for iBasin=1:671
    input(iBasin).precip=rainplusmelt(1:4018,iBasin);
    input(iBasin).pet=PET(1:4018,iBasin);
end


DPWMQ=zeros(4018,671);
DPWMET=zeros(4018,671);
for iBasin=1:671
    [fluxOutput,fluxInternal,storeInternal]=DPWM(input(iBasin),parameters(iBasin,:));
    DPWMQ(:,iBasin)=fluxOutput.Q;
    DPWMET(:,iBasin)=fluxOutput.ET;
end



