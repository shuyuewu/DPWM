


function [fluxOutput,fluxInternal,storeInternal]=DPWM(input,parameter)
% Copyright (C) 2024 Shuyue Wu
% Model reference
% Shuyue Wu,Jianshi Zhao, Murugesu Sivapalan. A Parsimonious Daily Water Balance Model based on the Proportionality Hypothesis.(submitted to Journal of Hydrology)

    P=input.precip;       % time series of rainfall and/or snowmelt          
    Ep=input.pet;         % time series of potential evapotranspiration
    lengthTS=length(P);   % length (i.e., num of days) of input time series 
    
    afa1=parameter(1);     % The ratio between initial catchment wetting and catchment wetting
    Smax=parameter(2);     % Soil moisture storage capacity
    afa2=parameter(3);     % The ratio between initial evapotranspiration opportunity and evapotranspiration opportunity
    kf=parameter(4);       % Fast runoff routing coefficient
    ks=parameter(5);       % Slow runoff routing coefficient, which should be smaller than kf
    
    Q=zeros(lengthTS,1);    % initialize time series of streamflow
    ET=zeros(lengthTS,1);   % initialize time series of evapotranspiration
    
    SM=zeros(lengthTS,1);   % initialize time series of soil moisture storage
    SS=zeros(lengthTS,1);   % initialize time series of fast routing storage
    SF=zeros(lengthTS,1);   % initialize time series of slow routing storage

    RF=zeros(lengthTS,1);  % initialize time series of fast runoff
    QF=zeros(lengthTS,1);  % initialize time series of fast streamflow
    RS=zeros(lengthTS,1);  % initialize time series of slow runoff
    QS=zeros(lengthTS,1);  % initialize time series of slow streamflow

    
    
    % Model warm-up periods are used to reduce the impact of uncertain initial conditions on model performance.
    warmSpan=365*1;  % defaultly use the first year for model warm-up  
    % determine the initial storages in an iterative procedure by 
    % letting the model repeat Year 1 of the data period until 
    % the stores reach an equilibrium for the first day of the year.
    SM00=0; 
    SF00=0;   
    SS00=0;
    dSM=10; 
    dSF=10;
    dSS=10; 
    it=0;
    
    while it<=1000&&(dSM>1 || dSS>1||dSF>1)%%
        it=it+1;
        SM0=SM00;
        SF0=SF00;
        SS0=SS00;
        
        for iDay=1:warmSpan
            A=P(iDay);                                                                    % available water for fisrt-stage of hydrological partitioning
            b=Smax-SM0+Ep(iDay);                                                          % catchment wetting capacity
            W=(A+b)/2/afa1/(2-afa1)-sqrt(((A+b)/2/afa1/(2-afa1))^2-A*b/afa1/(2-afa1));    % catchment wetting
            RF(iDay)=P(iDay)-W;                                                           % fast runoff
            QF(iDay)=kf/(1+kf)*(SF0+RF(iDay));                                            % fast streamflow
            SF(iDay,1)=SF0+RF(iDay)-QF(iDay);                                             % fast routing storage
            A=SM0+W;                                                                      % available for second-stage of hydrological partitioning
            b1=Ep(iDay)+Smax;                                                             % potential evapotranspiration opportunity
            Y=(A+b1)/2/afa2/(2-afa2)-sqrt(((A+b1)/2/afa2/(2-afa2))^2-A*b1/afa2/(2-afa2)); % evapotranspiration opportunity
            RS(iDay)=A-Y;                                                                 % slow runoff
            ET(iDay,1)=Y*Ep(iDay)/b1;                                                     % evapotranspiration
            SM(iDay,1)=Y-ET(iDay);                                                        % soil moisture storage at the end of the day
            QS(iDay)=ks/(1+ks)*(SS0+RS(iDay));                                            % slow streamflow
            SS(iDay,1)=SS0+RS(iDay)-QS(iDay);                                             % slow routing storage
            Q(iDay,1)=QF(iDay)+QS(iDay);                                                  % streamflow
            SM0=SM(iDay);
            SF0=SF(iDay);
            SS0=SS(iDay); 
        end
        dSM=abs(SM(warmSpan)-SM00);
        dSF=abs(SF(warmSpan)-SF00);
        dSS=abs(SS(warmSpan)-SS00);
        
        SM00=SM(warmSpan);
        SF00=SF(warmSpan);
        SS00=SS(warmSpan);
    end
    
    % water balance modeling for days after the warm-up period
    SM0=SM00;  
    SF0=SF00;
    SS0=SS00;  
    for iDay=warmSpan+1:lengthTS
        A=P(iDay);
        b=Smax-SM0+Ep(iDay);
        W=(A+b)/2/afa1/(2-afa1)-sqrt(((A+b)/2/afa1/(2-afa1))^2-A*b/afa1/(2-afa1));
        RF(iDay)=P(iDay)-W;%%
        QF(iDay)=kf/(1+kf)*(SF0+RF(iDay));%%
        SF(iDay,1)=SF0+RF(iDay)-QF(iDay);%% 
        A=SM0+W;
        b1=Ep(iDay)+Smax;
        Y=(A+b1)/2/afa2/(2-afa2)-sqrt(((A+b1)/2/afa2/(2-afa2))^2-A*b1/afa2/(2-afa2));
        RS(iDay)=A-Y;
        ET(iDay,1)=Y*Ep(iDay)/b1;
        SM(iDay,1)=Y-ET(iDay);
        QS(iDay)=ks/(1+ks)*(SS0+RS(iDay));
        SS(iDay,1)=SS0+RS(iDay)-QS(iDay);
        Q(iDay,1)=QF(iDay)+QS(iDay);
        SM0=SM(iDay);
        SF0=SF(iDay);%
        SS0=SS(iDay);        
        
    end
    fluxOutput.ET=ET;
    fluxOutput.Q=Q;
    storeInternal.SM=SM;
    storeInternal.SF=SF;
    storeInternal.SS=SS;
    fluxInternal.RF=RF;
    fluxInternal.QF=QF;
    fluxInternal.RS=RS;
    fluxInternal.QS=QS;
    
 
    
    
end
