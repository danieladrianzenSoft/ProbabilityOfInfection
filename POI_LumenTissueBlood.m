
function POI_LumenTissueBlood(params, makePlots)

Vir = Virus;
Vir.diffCoeffLumen = 1.27*10^(-8); %D_G
Vir.diffCoeffEpithelium = 1*10^(-10); %D_E previously 9x10^(-10)
Vir.diffCoeffStroma = 2*10^(-9); %D_S %previously 2x10^(-9)
Vir.lossDilutionInLumen = 0.27/3600; %k_D
Vir.lossStromaToBlood = params.k_B; %k_B -------------------------------------------------- %%% THIS IS COMPLETELY UNKNOWN!!!!
Vir.lossClearanceInBlood = params.k_B; %k_L
Vir.partitionLumenEpthelium = 0.3; %phi_GE %Hope et al, 2015. In Macaques, and Humans
Vir.partitionEpitheliumStroma = 1; %phi_ES Estimated
Vir.volumeOfDistribution = 5*1000; %volume of blood in ml
Vir.initialLumenConcentration = params.V_0; %V_0 in virions/ml
Vir.radius = 100*10^(-7)/2; %in cm
Vir.ratioCollisionsCausingInfection = 1/10; % 10%
Vir.infectivity = (6.5*10^(-7))/(24*3600); %perelson paper

Dru = Drug;
Dru.diffCoeffLumen = 6*10^(-6); %D_G cm2/s
Dru.diffCoeffEpithelium = 7*10^(-8); %D_E cm2/s
Dru.diffCoeffStroma = 4*10^(-7); %D_S cm2/s
Dru.lossDilutionInLumen = 1.22/3600; %k_D
Dru.lossStromaToBlood = 0.119/3600; %k_B
Dru.lossClearanceInBlood = 1.41/3600; %k_L
Dru.partitionLumenEpthelium = 0.75; %phi_GE
Dru.partitionEpitheliumStroma = 1; %phi_ES
Dru.volumeOfDistribution = 75*1000; %Vb ml
Dru.initialLumenConcentration = params.C_G0; %C_0 ng/ml
Dru.rateActivation = log(2)/3600; %kon
Dru.rateDeactivation = log(2)/(7*24*3600); %Koff
Dru.ratioActivationFromDeactivated = 0.1;

Tar = TCell;
Tar.initialStromaConcentration = 5*10^4; %T_0 Immunological Microenvironments in the Human Vagina, pudney, anderson 2005. ------ CHECK THIS!
%Tar.initialBloodConcentration = 1e6; %T_0 PBMC
Tar.initialBloodConcentration = 1e6; %T_0 PBMC
Tar.volumeFractionEpithelium = 0.95; %phiDP_E = 0.95; jing 2D
Tar.volumeFractionStroma = 0.1; %phiDP_S = 0.1; jing 2D 
Tar.deathRateStroma = 0.01/(24*3600); %s^(-1)
Tar.deathRateBlood = 0.01/(24*3600); %s^(-1)
Tar.productionRateStroma = Tar.deathRateStroma * Tar.initialStromaConcentration;
Tar.productionRateBlood = Tar.deathRateBlood * Tar.initialBloodConcentration;
Tar.radius = 10*10^(-4)/2; %diameter = 10um

Inc = ICell;
Inc.productionRateVirus = params.rho; %rho
Inc.deathRateStroma = 0.39/(24*3600); %dI stroma
Inc.deathRateBlood = 0.39/(24*3600); %dI blood

Lat = LCell;
Lat.deathRateStroma = 0.001/(24*3600); %dL stroma %Libing Rong, Perelson, 2009
Lat.deathRateBlood = 0.001/(24*3600); %dL blood
Lat.ratioInfectionsResultingInLatency = 0.001; %eta
Lat.transitionRateLtoI = 0.1/(24*3600); %a
Lat.proliferationRate = 0.015/(24*3600); %lambda_L

Geometry = Vagina;
Geometry.thicknessLumen = 0.04; %h_G cm
Geometry.thicknessEpithelium = params.h_E; %h_E cm
Geometry.thicknessStroma = 0.28; %h_S cm
Geometry.width = 3.35; % w cm
Geometry.length = 13; %l cm

Mesh = Discretization1D;
Mesh.numElements = 800;
Mesh.numLumen = ceil((Geometry.thicknessLumen/(Geometry.thicknessLumen + ...
    Geometry.thicknessEpithelium + Geometry.thicknessStroma)) * Mesh.numElements);
Mesh.numEpithelium = ceil((Geometry.thicknessEpithelium/(Geometry.thicknessLumen + ...
    Geometry.thicknessEpithelium + Geometry.thicknessStroma)) * Mesh.numElements);
Mesh.numStroma = Mesh.numElements - (Mesh.numLumen + Mesh.numEpithelium);
Mesh.x = [linspace(0, Geometry.thicknessLumen, Mesh.numLumen), ...
          linspace(Geometry.thicknessLumen+(Geometry.thicknessEpithelium/Mesh.numEpithelium), ...
            Geometry.thicknessLumen + Geometry.thicknessEpithelium, Mesh.numEpithelium) ...
          linspace(Geometry.thicknessLumen + Geometry.thicknessEpithelium + (Geometry.thicknessStroma / Mesh.numStroma), ...
            Geometry.thicknessLumen + Geometry.thicknessEpithelium + Geometry.thicknessStroma, Mesh.numStroma)];

% realx = [linspace(0,h_G,indG), linspace(h_G+(h_E/indE), h_G+h_E,indE) ...
%          linspace(h_G+h_E+(h_S/indS), h_G+h_E+h_S, indS)];
% T_0 = 10^(4)*1.088*10^(3)*1/100; % *1.088*10^(3)*1/100 cells/ml
% T_0blood = 2e6; %cells/ml
% %T_0 = 10^(4)*1.088*10^(3)*1/100; %cells/mg of tissue converted to cells/ml, given density of 1.088g/ml, and given 1% are cd4+ (perelson seminal)
% I_0 = 0;
% L_0 = 0;
% V_0 = params.V_0;
% 
% %delay = -2 * (60 * 60);
delay = params.T_VD;
% %delay = 0;
% %C_G0 = (1 * 10^(7));
% C_G0 = params.C_G0;
% %V_0 = (1 * 10^(4)); %Virions / ml, Katz lectures XXXXXXXXXX
% %V_0 = (1 * 10^(4));
% %filename = 'POI_LongTerm.csv';
% 
% %%%%%% PARAMETERS %%%%%%
% 
% %%%VIRAL DYNAMICS
% dt = 0.01/(24*3600); %s^(-1) death rate of target cells XXXXXXXXXX
% %dt = 0.02/(24*3600);
% lambda = dt*T_0;
% %lambda = 10^(4)/(24*3600); %cells/mL/s production rate of target cells
% %lambda = 100/(24*3600); %cells/ml/day
% %rho = 1.5*10^(3)/(24*3600); %s^(-1) production rate of V from I
% %rho = 850/(24*3600); 
% rho = params.rho;
% %rho = 1800/(24*3600); %XXXXXXXXX
% %rho = 1*10^(3)/(24*3600); %s^(-1) production rate of V from I
% %rho = 0;
% %del = 0.5/(24*3600); %s^(-1) death rate of I.
% %del = 0.80/(24*3600);
% del = 0.39/(24*3600);  %s^(-1) death rate of I. XXXXXXXXXXXX
% %beta = params.beta; 
% betaInf = (6.5*10^(-7))/(24*3600); %XXXXXXXXXXXXX
% %w = 10^(-3)/(24*3600); %s^(-1) per infected cell. Viral transmission
% %(cell-cell) XXXXXXXXXX
% w = 0;
% 
% %%%Latent cells:
% eta = 0.1; %fraction of infections that result in latency
% %eta = 10^(-3);
% %eta = 0;
% %eta = 0;
% dl = 0.001/(24*3600); %death rate of latently infected cells
% lambdal = 0.0045/(24*3600); %%NOT BEING USED RIGHT NOW 
% %a = 0.05/(24*3600); %transition rate of latent cells to productive infected cells.
% %a = 10^(-3)/(24*3600);
% a = 0.1/(24*3600);
% %a = 0;
% 
% %%%VIRAL TRANSPORT
% %D_G = 1.27*10^(-8); %XXXXXXXX
% D_G = 1.27*10^(-8);
% %D_G = 6*10^(-6);
% %D_E = 3*10^(-10); %previously 9x10^(-10)XXXXXXXXX
% D_E = 1*10^(-10);
% %D_E = 9*10^(-10);
% %D_S = 2*10^(-9); %previously 2x10^(-9)XXXXXXXXXX
% D_S = 1*10^(-9);
% %D_S = params.D_vS;
% %D_S = 4*10^(-7);
% k_D = 0.27/3600;
% %k_D = 1.22/3600; %XXXXXXXXXX
% %k_D = 0.551/3600;
% %k_b = 3/(24*3600);
% k_b = params.k_B;
% %k_b = 10/(24*3600); %previously 15
% %k_b = 15/(24*3600); %XXXXXXXXXX
% c = k_b; %s^(-1) elimination/clearance in tissue
% %Partition coefficient V_T/V_G
% %phiGE = 1;
% 
% kv_L = c;
% %phiGE = ((1.49*10^(4))/(4*10^(6))); %Carias, Hope et al. 2013;
% phiGE = 0.3; %Hope et al, 2015. In Macaques, and Humans
% phiES = 1;
% Vvb = 5*1000; %volume of blood, here roughly 5L, converted to ml.
% 
% %%%TFV TRANSPORT
% Dd_G = 6*10^(-6);
% %D_G = 3.72*10^(-6);
% %D_E = 5.73*10^(-8);
% %D_E = 5*10^(-8);
% Dd_E = 7*10^(-8); %PREVIOUSLY 7
% %D_S = 5.73*10^(-9);
% %D_S = 3.31*10^(-7); %Diffusion coefficient in stroma, from Chapter 5.
% %D_S = 5*10^(-8);
% Dd_S = 4*10^(-7);%PREVIOUSLY 
% phid_GE = 0.75; 
% phid_ES = 1;
% %k_b = 0;
% kd_D = 1.22/3600;
% kd_B = 0.119/3600;
% %kd_B = 100/3600;
% kd_L = 1.41/3600;
% Vdb = 75*1000;
%  
% %%%TFV-DP PRODUCTION:
% Kon = log(2)/3600;
% Koff = log(2)/(7*24*3600);
% r = 0.1;
% %These are volume fractions, not partition coefficients:
% phiDP_E = 0.95;
% phiDP_S = 0.1;
% 
% %Spatial and time constraints
% h_G = 0.04; %cm
% %h_E = 0.02; %cm
% h_E = params.h_E;
% h_S = 0.28; %cm
% 
% W = 3.35; %cm
% L = 13; %cm
% %V = 4;
% 
% numx = 800;
% 
% indG  = ceil((h_G/(h_G+h_E+h_S))*numx);
% indE  = ceil(((h_G + h_E) / (h_G + h_E + h_S)) * numx) - ...
%        (ceil((h_G / (h_G + h_E + h_S) * numx)));
% indS  = numx-(indG+indE);
% 
% realx = [linspace(0,h_G,indG), linspace(h_G+(h_E/indE), h_G+h_E,indE) ...
%          linspace(h_G+h_E+(h_S/indS), h_G+h_E+h_S, indS)];

% realx = linspace(0,h_G+h_E+h_S,numx);
% indG  = find(realx>h_G,1)
% indE  = find(realx>h_E,1)
% indS  = numx-(indG+indE)



%Virus: Gel, Epi, Stro, 
%Target cells: Stro, 
%Latent cells: Stro,
%Infected cells: Stro
%TFV: Gel, Epi, Stro
%TFV-DP: Epi, Stro

%totLength = h_G+h_E+h_S;

% x = [realx,totLength+realx(2:end),2*totLength+realx(indG+1:indG+indE+indS)-realx(indG)];

%Virus, TFV, TFV-DP, T, L, I

%x = [realx,realx,realx(indG+1:indE+indS),realx(indG+indE+1:indG+indE+indS),realx(indG+indE+1:indG+indE+indS),realx(indG+indE+1:indG+indE+indS)];

%x = [realx,totLength+realx(indG+indE+1:indG+indE+indS)-realx(indG+indE),totLength+h_S+realx(indG+indE+1:indG+indE+indS)-realx(indG+indE),...
%   totLength+2*h_S+realx(indG+indE+1:indG+indE+indS)-realx(indG+indE),totLength+3*h_S+realx,...
%   2*totLength+3*h_S+realx(indG+1:indG+indE+indS)-realx(indG)];
S = calcsparsity(Mesh.numElements,Mesh.numLumen,Mesh.numEpithelium,Mesh.numStroma);
S = sparse(S);

%opts1 = odeset('JPattern',S);
opts1 = odeset('Vectorized','on','JPattern',S);
finalT = 60*24*60*60; %100 days
numT = 1000;
%tau_0 = 13*24*60*60; %intracellular delay
tau_0 = 0;

%[f,maxPot,minPot,maxConc,minConc] = processingDRCurve(T_0);

if (Dru.initialLumenConcentration == 0)
    
    t = linspace(0,finalT,numT);
    tspan = [min(t),max(t)];
    
    IC = [Vir.initialLumenConcentration.*ones(1,Mesh.numLumen) 0.*ones(1,Mesh.numEpithelium) 0.*ones(1,Mesh.numStroma)... %virions
        Dru.initialLumenConcentration.*ones(1,Mesh.numLumen) 0.*ones(1,Mesh.numEpithelium) 0.*ones(1,Mesh.numStroma)... %TFV ng/ml
        0.*ones(1,Mesh.numEpithelium) 0.*ones(1,Mesh.numStroma),... %TFV-DP
        Tar.initialStromaConcentration.*ones(1,Mesh.numStroma) 0.*ones(1,Mesh.numStroma) 0.*ones(1,Mesh.numStroma),... %cells (target, latent, infected)
        0, 0, Tar.initialBloodConcentration, 0, 0, 0]; %TFV, TFV-DP, T, L, I, V in blood
    %ICtemp = [1*ones(1,indG) 1.*ones(1,indE) 1.*ones(1,indS)... %virions
    %    0.*ones(1,indG) 0.*ones(1,indE) 0.*ones(1,indS)... %TFV ng/ml
    %    1.*ones(1,indE) 1.*ones(1,indS),... %TFV-DP
    %    0.*ones(1,indS) 1.*ones(1,indS) 0.*ones(1,indS)]; %cells (target, latent, infected)

    
    %tic
    [t,V] = ode15s(@(t,V) dvdt(t,V,IC,Mesh,Geometry,... %Geometry
            Vir,... %Virus
            Tar, Inc, Lat,... %Cells
            Dru,... %TFV Drug
            tau_0), tspan, IC, opts1); %Intracellular delay
    %toc
    
    V = V';
    
elseif (delay == 0 && Dru.initialLumenConcentration ~=0)


    Geometry.thicknessLumen = 2*Geometry.thicknessLumen;
    Dru.initialLumenConcentration = Dru.initialLumenConcentration/2;
    Vir.initialLumenConcentration = Vir.initialLumenConcentration/2;
    %h_G1 = 2*Geometry.thicknessLumen;
    %C_G01 = D.initialLumenConcentration/2;
    %V_01 = V.initialLumenConcentration/2;
%     h_G1 = h_G;
%     C_G01 = C_G0;
%     V_01 = V_0;
    
    IC = [Vir.initialLumenConcentration.*ones(1,Mesh.numLumen) 0.*ones(1,Mesh.numEpithelium) 0.*ones(1,Mesh.numStroma)... %virions
        Dru.initialLumenConcentration.*ones(1,Mesh.numLumen) 0.*ones(1,Mesh.numEpithelium) 0.*ones(1,Mesh.numStroma)... %TFV ng/ml
        0.*ones(1,Mesh.numEpithelium) 0.*ones(1,Mesh.numStroma),... %TFV-DP
        Tar.initialStromaConcentration.*ones(1,Mesh.numStroma), 0.*ones(1,Mesh.numStroma), 0.*ones(1,Mesh.numStroma)... %cells (target, latent, infected)
        0, 0, Tar.initialBloodConcentration, 0, 0, 0]; %TFV, TFV-DP, T, L, I, V in blood
    t = linspace(delay,finalT,numT);
    tspan = [min(t),max(t)];
    
    %tic
    [t,V] = ode15s(@(t,V) dvdt(t,V,IC,Mesh,Geometry,... %Geometry
            Vir,... %Virus
            Tar, Inc, Lat,... %Cells
            Dru,... %TFV Drug
            tau_0), tspan, IC, opts1); %Intracellular delay
    %toc
    
    V = V';

elseif (delay < 0 && Dru.initialLumenConcentration ~=0)
    
    h_G1 = h_G; 
    C_G01 = C_G0;
    
    t1 = linspace(0,-1*delay,numT);
    tspan1 = [min(t1),max(t1)];

    IC1 = [0.*ones(1,Mesh.numLumen) 0.*ones(1,Mesh.numEpithelium) 0.*ones(1,Mesh.numStroma)... %virions
        Dru.initialLumenConcentration.*ones(1,Mesh.numLumen) 0.*ones(1,Mesh.numEpithelium) 0.*ones(1,Mesh.numStroma)... %TFV ng/ml
        0.*ones(1,Mesh.numEpithelium) 0.*ones(1,Mesh.numStroma),... %TFV-DP
        Tar.initialStromaConcentration.*ones(1,Mesh.numStroma), 0.*ones(1,Mesh.numStroma), 0.*ones(1,Mesh.numStroma) ... %cells (target, latent, infected)
        0, 0, Tar.initialBloodConcentration, 0, 0, 0]; %TFV, TFV-DP, T, L, I, V in blood
    
    [t1,V1] = ode15s(@(t1,V1) dvdt(t1,V1,IC1,Mesh,Geometry,... %Geometry
            Vir,... %Virus
            Tar, Inc, Lat,... %Cells
            Dru,... %TFV Drug
            tau_0), tspan1, IC1, opts1); %Intracellular delay
    
    %Later, semen deposition, so h_G to 0.08cm, and 1:1 dilution -> C_G/2, V_0/2
    Vir.initialLumenConcentration = Vir.initialLumenConcentration/2;
    Geometry.thicknessLumen = 2*Geometry.thicknessLumen;
    %V_01 = V_0/2;
    %h_G1 = 2*h_G; 

    IC = [Vir.initialLumenConcentration.*ones(1,Mesh.numLumen) 0.*ones(1,Mesh.numEpithelium) 0.*ones(1,Mesh.numStroma)... %virions
        V1(end,Mesh.numElements+1:Mesh.numElements+Mesh.numLumen)/2 V1(end,Mesh.numElements+Mesh.numLumen+1:Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma)... %TFV ng/ml
        V1(end,Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma+1:Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma+Mesh.numEpithelium+Mesh.numStroma),... %TFV-DP
        Tar.initialStromaConcentration.*ones(1,Mesh.numStroma), 0.*ones(1,Mesh.numStroma), 0.*ones(1,Mesh.numStroma), ... %cells (target, latent, infected) 
        0, 0, Tar.initialBloodConcentration, 0, 0, 0]; %TFV, TFV-DP, T, L, I, V in blood
    
    %t2 = linspace(-1*delay,finalT,numT/2);
    t2 = linspace(0,finalT,numT);

    tspan2 = [min(t2),max(t2)];

    %tic
    [t2,V] = ode15s(@(t2,V) dvdt(t2,V,IC,Mesh,Geometry,... %Geometry
            Vir,... %Virus
            Tar, Inc, Lat,... %Cells
            Dru,... %TFV Drug
            tau_0), tspan2, IC, opts1); %Intracellular delay
    %toc
    
    %V = [V1;V]';
    %t = [t1;t2];
    V = V';
    t = t2;

elseif (delay > 0 && Dru.initialLumenConcentration ~= 0)
    
    h_G1 = h_G; 
    V_01 = V_0;
    
    ratio = 1/10;
    
    t1 = linspace(0,delay,numT*ratio);
    %numT*ratio
    %(1-ratio)*numT
    tspan1 = [min(t1),max(t1)];
    
    IC1 = [Vir.initialLumenConcentration.*ones(1,Mesh.numLumen) 0.*ones(1,Mesh.numEpithelium) 0.*ones(1,Mesh.numStroma)... %virions
        0.*ones(1,Mesh.numLumen) 0.*ones(1,Mesh.numEpithelium) 0.*ones(1,Mesh.numStroma)... %TFV ng/ml
        0.*ones(1,Mesh.numEpithelium) 0.*ones(1,Mesh.numStroma),... %TFV-DP
        Tar.initialStromaConcentration.*ones(1,Mesh.numStroma), 0.*ones(1,Mesh.numStroma), 0.*ones(1,Mesh.numStroma)... %cells (target, latent, infected)
        0, 0, Tar.initialBloodConcentration, 0, 0, 0]; %TFV, TFV-DP, T, L, I, V in blood

    [t1,V1] = ode15s(@(t1,V1) dvdt(t1,V1,IC1,Mesh,Geometry,... %Geometry
            Vir,... %Virus
            Tar, Inc, Lat,... %Cells
            Dru,... %TFV Drug
            tau_0), tspan1, IC1, opts1); %Intracellular delay
    
    Dru.initialLumenConcentration = Dru.initialLumenConcentration / 2;
    Geometry.thicknessLumen = Geometry.thicknessLumen * 2;
    %C_G01 = C_G0/2;
    %h_G1 = h_G*2; 
%     h_G1 = h_G;
%     C_G01 = C_G0;
    
    IC = [V1(end,1:Mesh.numLumen)/2 V1(end,Mesh.numLumen+1:Mesh.numElements)... %virions
        Dru.initialLumenConcentration.*ones(1,Mesh.numLumen) 0.*ones(1,Mesh.numEpithelium) 0.*ones(1,Mesh.numStroma)... %TFV ng/ml
        0.*ones(1,Mesh.numEpithelium) 0.*ones(1,Mesh.numStroma),... %TFV-DP
        Tar.initialStromaConcentration.*ones(1,Mesh.numStroma), 0.*ones(1,Mesh.numStroma), 0.*ones(1,Mesh.numStroma)... %cells (target, latent, infected)
        0, 0, T.initialBloodConcentration, 0, 0, 0]; %TFV, TFV-DP, T, L, I, V in blood
    
    t2 = linspace(delay,finalT,(1-ratio)*numT);
    tspan2 = [min(t2),max(t2)];

    %tic
    [t2,V] = ode15s(@(t2,V) dvdt(t2,V,IC,Mesh,Geometry,... %Geometry
            Vir,... %Virus
            Tar, Inc, Lat,... %Cells
            Dru,... %TFV Drug
            tau_0), tspan2, IC, opts1); %Intracellular delay
    %toc
    
    V = [V1;V]';
    t = [t1;t2];
    
end

% %CALCULATING MICROBICIDE EFFICACY Q
% 
% %Original data from 
% TFVDPPotency = load('TFVDP_Potency.txt');
% smoothTFVDPPotency = smooth(TFVDPPotency(:,1),TFVDPPotency(:,2),200,'rloess');
% % ConcX = log10(10.^(TFVDPPotency(:,1))*(10^(6)/180)); %4.9917*10^3 10^6cells / mL tissue  = 4.9917 10^6cells / mg tissue. 
% %                                                %Assuming: density tissue = 1000mg/mL. Conversion from Schwartz et al.
% %                                                %fmol/10^6 cells to fmol/mg tissue
% ConcX = log10(10.^(TFVDPPotency(:,1))*(T_0/(10^(6))));
% %ConcX = log10(10.^(TFVDPPotency(:,1))*4.9917*10^(-2));
%                      
% f = fit(ConcX,smoothTFVDPPotency,'fourier8');
% fvals = f.a0 + f.a1*cos(ConcX*f.w) + f.b1*sin(ConcX*f.w) + ... %Fit for concentration in units of log10 fmol/mg tissue
%         f.a2*cos(2*ConcX*f.w) + f.b2*sin(2*ConcX*f.w) + ...
%         f.a3*cos(3*ConcX*f.w) + f.b3*sin(3*ConcX*f.w) + ...
%         f.a4*cos(4*ConcX*f.w) + f.b4*sin(4*ConcX*f.w) + ...
%         f.a5*cos(5*ConcX*f.w) + f.b5*sin(5*ConcX*f.w) + ...
%         f.a6*cos(6*ConcX*f.w) + f.b6*sin(6*ConcX*f.w) + ...
%         f.a7*cos(7*ConcX*f.w) + f.b7*sin(7*ConcX*f.w) + ...
%         f.a8*cos(8*ConcX*f.w) + f.b8*sin(8*ConcX*f.w);
% maxPotf = min(fvals);
% minPotf = max(fvals);
% maxConcf = max(ConcX);
% minConcf = min(ConcX);

%[qtemp,Pottemp] = PotCalc(ConcX,f,maxPotf,minPotf,maxConcf,minConcf);

%Drug
%CTFV = V(1:ng+ne+ns,:);
%CTFVDP = [zeros(ng,length(t2));C(ng+ne+ns+1:ng+ne+ne+ns+ns,:)];
Cg = V(Mesh.numElements:Mesh.numElements+Mesh.numLumen,:);
Ce = V(Mesh.numElements+Mesh.numLumen+1:Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium,:);
Cs = V(Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+1:Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma,:);
Cdpe = V(Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma+1:Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma+Mesh.numEpithelium,:);
Cdps = V(Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma+Mesh.numEpithelium+1:Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma+Mesh.numEpithelium+Mesh.numStroma,:);
Cb = V((Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma)+Mesh.numEpithelium+4*Mesh.numStroma+1,:);
Cbdp = V((Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma)+Mesh.numEpithelium+4*Mesh.numStroma+2,:);
Tb = V((Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma)+Mesh.numEpithelium+4*Mesh.numStroma+3,:);
Lb = V((Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma)+Mesh.numEpithelium+4*Mesh.numStroma+4,:);
Ib = V((Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma)+Mesh.numEpithelium+4*Mesh.numStroma+5,:);
Vblood = V((Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma)+Mesh.numEpithelium+4*Mesh.numStroma+5,:);

Cdps_FM = (Cdps/447.173)*10^(3); %fmol/mg
%Ct = C(ng+1:ng+ne+ns,:);
CG_avg = trapz(Cg)/(Mesh.numLumen-1);
CE_avg = trapz(Ce)/(Mesh.numEpithelium-1);
CS_avg = trapz(Cs)/(Mesh.numStroma-1);
CDP_Eavg = trapz(Cdpe)/(Mesh.numEpithelium-1);
CDP_Savg = trapz(Cdps)/(Mesh.numStroma-1);
% CDP_SavgFM = (CDP_Savg/447.173)*10^(3); %from ng/ml to fmol/mg
% CDP_Tavg = (CDP_Eavg*h_E+CDP_Savg*h_S)/(h_E+h_S);
% CT_avg = (CE_avg*h_E+CS_avg*h_S)/(h_E+h_S);

%Virions & Cells
vir = V(1:Mesh.numElements,:);
tar = V((Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma)+Mesh.numEpithelium+Mesh.numStroma+1:(Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma)+Mesh.numEpithelium+2*Mesh.numStroma,:);
%size(q)
%size(V(numx+1:numx+indS,:))
%tarDr = (1-q).*V((numx+indG+indE+indS)+indE+indS+1:(numx+indG+indE+indS)+indE+2*indS,:);
lat = V((Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma)+Mesh.numEpithelium+2*Mesh.numStroma+1:(Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma)+Mesh.numEpithelium+3*Mesh.numStroma,:);
inf = V((Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma)+Mesh.numEpithelium+3*Mesh.numStroma+1:(Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma)+Mesh.numEpithelium+4*Mesh.numStroma,:); 
Vg = V(1:Mesh.numLumen,:);
Ve = V(Mesh.numLumen+1:Mesh.numLumen+Mesh.numEpithelium,:);
Vs = V(Mesh.numLumen+Mesh.numEpithelium+1:Mesh.numElements,:);
Vt = V(Mesh.numLumen+1:Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma,:);

TT_avg = (trapz(tar,1)/(Mesh.numStroma-1))';
IT_avg = (trapz(inf,1)/(Mesh.numStroma-1))';
LT_avg = (trapz(lat,1)/(Mesh.numStroma-1))';

% collisionsVC_t= 0.03*4*pi*W*L*(Rc+Rv)*D_S*collisionsVCIntX;
% collisionsVC_cum = 0.03*4*pi*W*L*(Rc+Rv)*D_S*cumtrapz(t,collisionsVCIntX);
% collisionsMC_t = 4*pi*W*L*(Rc+Rm)*Dd_S*collisionsMCIntX;
% collisionsMC_cum = 4*pi*W*L*(Rc+Rm)*Dd_S*cumtrapz(t,collisionsMCIntX);

%collisionsTest2 = cumtrapz(t,collisionsX);
%A = 4*pi*W*L*(Rc+Rv)*D_S
%At = 4*pi*(Rc+Rv)*D_S
VG_avg = (trapz(Vg,1)/(Mesh.numLumen-1))';
VE_avg = (trapz(Ve,1)/(Mesh.numEpithelium-1))';
VS_avg = (trapz(Vs,1)/(Mesh.numStroma-1))';
VT_avg = (VE_avg*Geometry.thicknessEpithelium+VS_avg*Geometry.thicknessStroma)/(Geometry.thicknessEpithelium+Geometry.thicknessStroma);


%B0 = [T_0blood;0;0;0;0;0];
%B0 = [0; 0];
% B0 = 0;
% tspan = [min(t) max(t)];
% opts1 = odeset('Vectorized','on');
% 
% params = { W, L, h_S, c, lambda, beta, dt, eta, dl, a, del, Vvb, kv_L, rho, kd_B, kd_L, Vdb, Kon, Koff, r, T_0blood };

%size(t)
%size(CS_avg)

% VS_vec = [t, VS_avg];
% CS_vec = [t, CS_avg'];

% figure()
% plot(t/3600,Cb,'LineWidth',2);
% xlim([0,24])

%[t,B] = ode15s(@(t,B) dbcdt(t,B,VS_vec,CS_vec,params),t ,B0);

%result = dbcdt(t,B,VS_vec,CS_vec,params);

%size(result)

if makePlots == 1
    realxTissue = Mesh.x(Mesh.numLumen+1:Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma)-0.04;
    TFVdpTissue = V(Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma+1:Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma+Mesh.numEpithelium+Mesh.numStroma,:);
    TFVTissue = V(Mesh.numElements+Mesh.numLumen+1:Mesh.numElements+Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma,:);
    TFVdp_avgs_Tissue = [CDP_Eavg', CDP_Savg'];
    TFV_avgs_Tissue = [CG_avg', CE_avg', CS_avg'];
    HIVTissue = V(Mesh.numLumen+1:Mesh.numLumen+Mesh.numEpithelium+Mesh.numStroma,:)+10^(-12);
    HIV_avgs_Tissue = [VG_avg, VE_avg, VS_avg, VT_avg];
    BloodComponents = [Cb',Cbdp',Tb', Lb', Ib', Vblood'];
    cellsTissue = [TT_avg,LT_avg,IT_avg];
    createLinePlots(t/3600,realxTissue,TFV_avgs_Tissue,TFVdp_avgs_Tissue,HIV_avgs_Tissue,cellsTissue,Mesh.numEpithelium,Mesh.numStroma,Geometry.thicknessEpithelium,Geometry.thicknessStroma);
    createHeatmap(t/3600,realxTissue,TFVdpTissue,HIVTissue,Geometry.thicknessEpithelium,Geometry.thicknessStroma)
    createBloodPlots(t/3600, BloodComponents);
end

% figure()
% plot(t/(24*3600),B(:,1), 'color', 'k', 'LineWidth',2)
% xlabel('Time (days)','FontSize',16,'FontWeight','Bold')
% ylabel('Target Cells','FontSize',16,'FontWeight','Bold')
% figure()
% plot(t/(24*3600),B(:,2), 'color', 'r', 'LineWidth',2)
% xlabel('Time (days)','FontSize',16,'FontWeight','Bold')
% ylabel('Latent Cells','FontSize',16,'FontWeight','Bold')
% figure()
% plot(t/(24*3600),B(:,3), 'color', 'b', 'LineWidth',2)
% xlabel('Time (days)','FontSize',16,'FontWeight','Bold')
% ylabel('Infected Cells','FontSize',16,'FontWeight','Bold')
% figure()
% plot(t/(24*3600),B(:,4), 'color', 'g', 'LineWidth',2)
% xlabel('Time (days)','FontSize',16,'FontWeight','Bold')
% ylabel('Virus','FontSize',16,'FontWeight','Bold')



%[t,y] = ode45(@(t,y) vdp1(t,y),[0 20],[2; 0]);

%figure()
%plot(t,y)
% b=2;
% d=3;
% time_period = [0 10]; %  tmin,tmax
% initial = [0,10,0,0]; % initial cond [t,y]
% [t,y]= ode45(@(t,y)myode45function(t,y,b,d), time_period, initial); %function call
% figure(1)
% plot(t,y(:,1))
% axis([0 10 0 12])
% title('ODE solution: y(t)=2-2*exp(-0.2*t) ','Fontweight','bold','fontsize',14);


%     [t2,V] = ode15s(@(t2,V) dvdt(t2,V,IC,realx,h_G1,h_E,h_S,indG,indE,indS,numx,W,L,... %Geometry
%             D_G,D_E,D_S,k_D,k_b,phiGE,phiES,... %Viral Transport
%             lambda,dt,rho,c,del,w,eta,lambdal,dl,a,T_0,... %Viral Dynamics
%             Dd_G,Dd_E,Dd_S,phid_GE,phid_ES,kd_D,kd_B,kd_L,Vdb,... %TFV Transport
%             Kon, Koff,r, phiDP_E,phiDP_S,tau_0), tspan2, IC,opts1);
%grad = gradient(VT_avg,t);

%tm20 = find(t>=(t(end)-(20*24*3600)),1);
% infected = trapz(tar,1)<trapz(T_0*ones(indS),1);

%% BLOOD COMPARTMENT DRUG

% cb_TFV = zeros(1,length(t));
% 
% f = CS_avg*W*L*h_S*kd_B*2;
% max(f)
% decay = exp(-t*kd_L);
% 
% for i = 2:length(t)
%     cb_TFV(i) = trapz(t(1:i),f(1:i).*decay(i:-1:1)/Vdb);
% end


%%

% infected = VS_avg(end)>0.01*max(VS_avg);
% infectionTP = find((TT_avg)<T_0,1);
% 
% if infected
%     t_inf = t(infectionTP);
% else
%     t_inf = 0;
% end

%     if isempty(infectionTP)
%         infected = 0;
%         t_inf = 0;
%     else
%         infected = 1;
%         t_inf = t(infectionTP);
%     end

end

function q = calcQ(cTFVDP)

IC50 = 180; %(this is in ng/ml, roughly = 0.4 uM, given molar mass of 447.173 g/mol))
%     if (IC(numx+1) ~= 0)
q = 1./(1+(IC50./(cTFVDP)));
        %q = q(:);
%     else
%         q = zeros(indS,1);
%     end

end


function DVDT = dvdt(t,V,IC,Mesh,Geometry,... %Geometry
            Vir,... %Virus
            Tar, Inc, Lat,... %Cells
            Dru,... %TFV Drug
            tau_0) %TFV-DP production

%NOTE: TO UNDERSTAND THESE EQUATIONS, LOOK UP FOR DIFFERENCE EQUATIONS: CENTRAL, BACKWARD
%AND FORWARD.

%m = 1;
%dxG = realx(2)-realx(1);
%dxE = realx(indG+2)-realx(indG+1);
%dxS = realx(indG+indE+2)-realx(indG+indE+1);

dVdt = zeros(length(IC),size(V,2));
%dVdx = zeros(length(IC),size(V,2));
d2Vdx2 = zeros(length(IC),size(V,2));

%sizeV = size(V)  

dxL = Geometry.thicknessLumen / (Mesh.numLumen - 1);
dxE = Geometry.thicknessEpithelium / (Mesh.numEpithelium - 1);
dxS = Geometry.thicknessStroma / (Mesh.numStroma - 1);

lossInEpithelium = 0/(24*3600);
beta_inf_epi = 0.1/100;
beta_inf_blood = (1e-8)/(24*3600);
decay_rate = 0.2/(24*3600);

%perCollisionsInfection = 1/10;
%Rc = 10*10^(-4)/2; %diameter = 10um
%Rv = 100*10^(-7)/2; %diameter = 100nm

Dv_L = Vir.diffCoeffLumen;
Dv_E = Vir.diffCoeffEpithelium;
Dv_S = Vir.diffCoeffStroma;
kv_D = Vir.lossDilutionInLumen;
kv_B = Vir.lossStromaToBlood;
kv_L = Vir.lossClearanceInBlood;
phiv_LE = Vir.partitionLumenEpthelium;
phiv_ES = Vir.partitionEpitheliumStroma;
Vvb = Vir.volumeOfDistribution;
Rv = Vir.radius;
s = Vir.ratioCollisionsCausingInfection;
betaInf = Vir.infectivity;
%s = beta_inf_epi + (Vir.ratioCollisionsCausingInfection - beta_inf_epi)*exp(-decay_rate*t);
%betaInf = beta_inf_blood + (Vir.infectivity - beta_inf_blood)*exp(-decay_rate*t);

Dd_L = Dru.diffCoeffLumen;
Dd_E = Dru.diffCoeffEpithelium;
Dd_S = Dru.diffCoeffStroma;
kd_D = Dru.lossDilutionInLumen;
kd_B = Dru.lossStromaToBlood;
kd_L = Dru.lossClearanceInBlood;
phid_LE = Dru.partitionLumenEpthelium;
phid_ES = Dru.partitionEpitheliumStroma;
Vdb = Dru.volumeOfDistribution;
Kon = Dru.rateActivation;
Koff = Dru.rateDeactivation;
r = Dru.ratioActivationFromDeactivated;

lambda_S = Tar.productionRateStroma;
lambda_B = Tar.productionRateBlood;
phiDP_E = Tar.volumeFractionEpithelium;
phiDP_S = Tar.volumeFractionStroma;
dT_S = Tar.deathRateStroma;
dT_B = Tar.deathRateBlood;
Rc = Tar.radius;
cellV = Tar.getVolume();

rho = Inc.productionRateVirus;
dI_S = Inc.deathRateStroma;
dI_B = Inc.deathRateBlood;

dL_S = Lat.deathRateStroma;
dL_B = Lat.deathRateBlood;
eta = Lat.ratioInfectionsResultingInLatency;
a = Lat.transitionRateLtoI;
lambda_L = Lat.proliferationRate;

%h_L = Geometry.thicknessLumen;
%h_E = Geometry.thicknessEpithelium;
h_S = Geometry.thicknessStroma;
W = Geometry.width;
L = Geometry.length;

numx = Mesh.numElements;
indL = Mesh.numLumen;
indE = Mesh.numEpithelium;
indS = Mesh.numStroma;


% dxG = h_G/(indG-1);
% dxE = h_E/(indE-1);
% dxS = h_S/(indS-1);

%HIV INTERFACE
V_intfLEa = (Dv_L.*V(indL,:)+(Dv_E.*V(indL+1,:)))./((Dv_E.*phiv_LE)+Dv_L); %C at gel/tissue interface
V_intfLEb = (Dv_L.*V(indL,:)+(Dv_E.*V(indL+1,:)))./((Dv_L./phiv_LE)+Dv_E);
V_intfESa = (Dv_E.*V(indL+indE,:)+(Dv_S.*V(indL+indE+1,:)))./((Dv_S.*phiv_ES)+Dv_E); %C at gel/tissue interface
V_intfESb = (Dv_E.*V(indL+indE,:)+(Dv_S.*V(indL+indE+1,:)))./((Dv_E./phiv_ES)+Dv_S);
% V_intf1a = (Vir.diffCoeffLumen.*V(Mesh.numLumen,:)+(Vir.diffCoeffEpithelium.*V(Mesh.numLumen+1,:)))./((D_E.*phiGE)+D_G); %C at gel/tissue interface
% V_intf1b = (D_G.*V(indG,:)+(D_E.*V(indG+1,:)))./((D_G./phiGE)+D_E);
% V_intf2a = (D_E.*V(indG+indE,:)+(D_S.*V(indG+indE+1,:)))./((D_S.*phiES)+D_E); %C at gel/tissue interface
% V_intf2b = (D_E.*V(indG+indE,:)+(D_S.*V(indG+indE+1,:)))./((D_E./phiES)+D_S);

%DRUG INTERFACE
C_intf_LEa = (Dd_L.*V(numx+indL,:)+(Dd_E.*V(numx+indL+1,:)))./((Dd_E.*phid_LE)+Dd_L); %Right before interface
C_intf_LEb = (Dd_L.*V(numx+indL,:)+(Dd_E.*V(numx+indL+1,:)))./((Dd_L./phid_LE)+Dd_E); %Right after interface
C_intf_ESa = (Dd_E.*V(numx+indL+indE,:)+(Dd_S.*V(numx+indL+indE+1,:)))./((Dd_S.*phid_ES)+Dd_E); %Right before interface
C_intf_ESb = (Dd_E.*V(numx+indL+indE,:)+(Dd_S.*V(numx+indL+indE+1,:)))./((Dd_E./phid_ES)+Dd_S); %Right after interface

%% VIRUS

ii = 1; %y direction BC zero flux (dcdx=0)
    dVdt(ii,:) = 2*Dv_L.*(V(ii+1,:)-V(ii,:))./(dxL.^2)-(kv_D.*V(ii,:));
    %dVdt(i) = 2*D_G.*(V(i+1)-V(i))./(dx.^2);
    
ii = 2:indL-1; %in gel
    %dVdx(ii,:) = (V(ii+1,:)-V(ii,:))./dxG;
    d2Vdx2(ii,:) = (V(ii+1,:)-2.*V(ii,:)+V(ii-1,:))./(dxL.^2);
    dVdt(ii,:) = Dv_L.*d2Vdx2(ii,:)-(kv_D.*V(ii,:));

ii = indL; %right before interface - gel
    %dVdx(ii,:) = (V_intf1a - V(ii,:))./dxG;
    d2Vdx2(ii,:) = (V(ii-1,:)-2.*V(ii,:)+V_intfLEa)./(dxL.^2);
    dVdt(ii,:) = Dv_L.*d2Vdx2(ii,:)-(kv_D.*V(ii,:));

ii = (indL+1); %right after interface -in epithelium
    %dVdx(ii,:) = (V(ii+1,:)-V(ii,:))./dxE;
    d2Vdx2(ii,:) = ((V_intfLEb)-2.*V(ii,:)+V(ii+1,:))./(dxE.^2);
    dVdt(ii,:) = Dv_E.*d2Vdx2(ii,:);

ii = (indL+2):(indL+indE)-1; %in epithelium
    %dVdx(ii,:) = (V(ii+1,:)-V(ii,:))./dxE;
    d2Vdx2(ii,:) = (V(ii+1,:)-2.*V(ii,:)+V(ii-1,:))./(dxE.^2);
    dVdt(ii,:) = Dv_E.*d2Vdx2(ii,:);

ii = indL+indE; %right before interface - epithelium,stroma
    %dVdx(ii,:) = (V_intf2a-V(ii,:))./dxE;
    d2Vdx2(ii,:) = (V(ii-1,:)-2.*V(ii,:)+V_intfESa)./(dxE.^2);
    dVdt(ii,:) = Dv_E.*d2Vdx2(ii,:);

ii = (indL+indE+1); %right after interface -in stroma
    %dVdx(ii,:) = (V(ii+1,:)-V(ii,:))./dxS;
    d2Vdx2(ii,:) = ((V_intfESb)-2.*V(ii,:)+V(ii+1,:))./(dxS.^2);
    dVdt(ii,:) = Dv_S.*d2Vdx2(ii,:) - kv_B.*V(ii,:) + rho*V(ii+(numx+indE+indS+indS+indS+indS),:);
    %dVdt(i) = D_S.*d2Vdx2(i) - k_b.*V(i);

ii = (indL+indE+2):numx-1; %in stroma
    %dVdx(ii,:) = (V(ii+1,:)-V(ii,:))./dxS;
    d2Vdx2(ii,:) = (V(ii+1,:)-2.*V(ii,:)+V(ii-1,:))./(dxS.^2);
    dVdt(ii,:) = Dv_S.*d2Vdx2(ii,:) - kv_B.*V(ii,:) + rho*V(ii+(numx+indE+indS+indS+indS+indS),:);
    %dVdt(i) = D_S.*d2Vdx2(i) - k_b.*V(i); 

ii = numx;%end of tissue
    %dCdt(i) = D_T.*(C(i)-2.*C(i-1))./(dx.^2) - k_b.*C(i); %BC zero conc
    %dCdt(i) = D_T.*(-2*C(i)+C(i-1))./(dx.^2) - k_b.*C(i); %BC zero conc CORRECT
    %dVdt(i) = 2*D_T*(V(i-1)-V(i))./(dx.^2) - k_b.*V(i); %BC zero flux
    %dVdt(:,i) = 2*D_S*(V(i-1,:)-V(i,:))./(dxS.^2) - k_b.*V(i,:) + rho*V((numx+indG+indE+indS)+indE+indS+3,:)/indS; %BC zero flux
    dVdt(ii,:) = 2*Dv_S*(V(ii-1,:)-V(ii,:))./(dxS.^2) - kv_B.*V(ii,:) + rho*V(ii+(numx+indE+indS+indS+indS+indS),:); %BC zero flux
    %dVdt(i) = 2*D_S*(V(i-1)-V(i))./(dx.^2) - k_b.*V(i); %BC zero flux

    %fprintf('%30s %.5f %s\n', 'Timing', toc(aaa2), 'seconds');
    %aaa3 = tic;

if (IC(numx+1) ~= 0) %%DRUG: SOLVE ONLY IF C_0 ~= 0:

    %% TFV
    ii = numx+1; %TFV at x = 0, zero flux B.C.
        %dVdt(i) = 2*Dd_G.*(V(i+1,:)-V(i,:))./(dxG.^2)-(kd_D.*V(i,:));
        %dVdt(i) = 2*Dd_G.*(V(i+1,:)-V(i,:))./(dxG.^2)-(kd_D.*V(i,:));
        dVdt(ii,:) = 2*Dd_L.*(V(ii+1,:)-V(ii,:))./(dxL.^2)-(kd_D.*V(ii,:));

    ii = numx+2:numx+indL-1; %TFV in gel
        %dVdx(ii,:) = (V(ii+1,:)-V(ii,:))./dxG;
        d2Vdx2(ii,:) = (V(ii+1,:)-2.*V(ii,:)+V(ii-1,:))./(dxL.^2);
        dVdt(ii,:) = Dd_L.*d2Vdx2(ii,:)-(kd_D.*V(ii,:));

    ii = numx+indL; %right before interface - gel/epithelium
        %dVdx(ii,:) = (C_intf_GEa - V(ii,:))./dxG;
        d2Vdx2(ii,:) = (V(ii-1,:)-2.*V(ii,:)+C_intf_LEa)./(dxL.^2);
        dVdt(ii,:) = Dd_L.*d2Vdx2(ii,:)-(kd_D.*V(ii,:));

    ii = (numx+indL+1); %right after interface - gel/epithelium
        %dVdx(ii,:) = (V(ii+1,:)-V(ii,:))./dxE;
        d2Vdx2(ii,:) = ((C_intf_LEb)-2.*V(ii,:)+V(ii+1,:))./(dxE.^2);
        dVdt(ii,:) = Dd_E.*d2Vdx2(ii,:);

    ii = (numx+indL+2):(numx+indL+indE-1); %in epithelium
        %dVdx(ii,:) = (V(ii+1,:)-V(ii,:))./dxE;
        d2Vdx2(ii,:) = (V(ii+1,:)-2.*V(ii,:)+V(ii-1,:))./(dxE.^2);
        dVdt(ii,:) = Dd_E.*d2Vdx2(ii,:);

    ii = numx+indL+indE; %right before interface - epithelium/stroma
        %dVdx(ii,:) = (C_intf_ESa - V(ii,:))./dxE;
        d2Vdx2(ii,:) = (V(ii-1,:)-2.*V(ii,:)+C_intf_ESa)./(dxE.^2);
        dVdt(ii,:) = Dd_E.*d2Vdx2(ii,:);

    ii = numx+indL+indE+1; %right after interface - epithelium/stroma
        %dVdx(ii,:) = (V(ii+1,:)-V(ii,:))./dxS;
        d2Vdx2(ii,:) = ((C_intf_ESb)-2.*V(ii,:)+V(ii+1,:))./(dxS.^2);
        dVdt(ii,:) = Dd_S.*d2Vdx2(ii,:)- kd_B.*V(ii,:);

    ii = numx+indL+indE+2:numx+indL+indE+indS-1; %in stroma
        %dVdx(ii,:) = (V(ii+1,:)-V(ii,:))./dxS;
        d2Vdx2(ii,:) = (V(ii+1,:)-2.*V(ii,:)+V(ii-1,:))./(dxS.^2);  
        dVdt(ii,:) = Dd_S.*d2Vdx2(ii,:) - kd_B.*V(ii,:); 

    ii = (numx+indL+indE+indS); %end of stroma
        %dCdt(i) = D_S.*(C(i)-2.*C(i-1))./(dx_S.^2) - k_B.*C(i); %BC zero conc
        %dCdt(i) = D_S.*(-2*C(i)+C(i-1))./(dx.^2) - k_B.*C(i); %BC zero conc CORRECT
        %dVdt(:,i) = 2*Dd_S*(V(i-1,:)-V(i,:))./(dxS.^2) - kd_B.*V(i,:); %BC zero flux
        dVdt(ii,:) = 2*Dd_S*(V(ii-1,:)-V(ii,:))./(dxS.^2) - kd_B.*V(ii,:); %BC zero flux
    %fprintf('%30s %.5f %s\n', 'Timing', toc(aaa3), 'seconds');
    %aaa4 = tic;

    %% TFV-DP

    ii = (numx+indL+indE+indS)+1:(numx+indL+indE+indS)+indE; %in epithelium
        bracks = (V(ii-indE-indS,:).*phiDP_E-V(ii,:)/r);
    %     if bracks < 0
    %         bracks = 0;
    %     end
        dVdt(ii,:) = Kon*bracks.*(bracks>0)-Koff*V(ii,:);

    ii = (numx+indL+indE+indS)+indE+1:(numx+indL+indE+indS)+indE+indS; %in stroma
        bracks = (V(ii-indE-indS,:).*phiDP_S-V(ii,:)/r);
    %     if bracks < 0
    %         bracks = 0;
    %     end
        dVdt(ii,:) = Kon*bracks.*(bracks>0)-Koff*V(ii,:);
        

%%%%TENOFOVIR CORRECTION DUE TO TENOFOVIR DIPHOSPHATE CONVERSION

    dVdt((numx+indL+1):(numx+indL+indE),:) = dVdt(numx+indL+1:numx+indL+indE,:) - dVdt((numx+indL+indE+indS)+1:(numx+indL+indE+indS)+indE,:); %in epithelium
    dVdt((numx+indL+indE+1):(numx+indL+indE+indS),:) = dVdt(numx+indL+indE+1:numx+indL+indE+indS,:) - dVdt((numx+indL+indE+indS)+indE+1:(numx+indL+indE+indS)+indE+indS,:); %in stroma

end 

%Vs = V(indG+indE+1:numx,:);
%tar = V((numx+indG+indE+indS)+indE+indS+1:(numx+indG+indE+indS)+indE+2*indS,:);
%inf = V((numx+indG+indE+indS)+indE+3*indS+1:(numx+indG+indE+indS)+indE+4*indS,:);

%collisionsCC_t = 4*pi*W*L*(Rc+Rc)*(10^-5)*tar*inf;

%aaa5 = tic;

MultConcCollisions = V((numx+indL+indE+indS)+indE+indS+1:(numx+indL+indE+indS)+indE+2*indS,:).*V((indL+indE+1):numx,:);
collisionsVC_t= 4*pi*(Rc+Rv)*(Dv_S)*MultConcCollisions;
cellInfections = s*collisionsVC_t;

%Cdps_uM = (V((numx+indG+indE+indS)+indE+1:(numx+indG+indE+indS)+indE+indS)/447.173);

% if (IC(numx+1) ~= 0)
%     %q = PotCalc(log10(Cdps_FM),f,maxPot,minPot,maxConc,minConc);
%     IC50 = 0.2; %in uM
%     q = 1./(1+(IC50./(Cdps_uM)));
% else
%     q = zeros(indS,1);
% end

% IC50 = 0.2; %uM
% 
% if (IC(numx+1) ~= 0)
%     q = 1./(1+(IC50./(V((numx+indG+indE+indS)+indE+1:(numx+indG+indE+indS)+indE+indS,:))));
%     %q = q(:);
% else
%     q = zeros(indS,1);
% end
cTFVDP = (V((numx+indL+indE+indS)+indE+1:(numx+indL+indE+indS)+indE+indS,:));

if (IC(numx+1) ~= 0)
    q = calcQ(cTFVDP);
else 
    q = zeros(indS,1);
end


% q = calcQ(cTFVDP, numx, indS, IC);

%fprintf('%30s %.5f %s\n', 'Timing', toc(aaa5), 'seconds');

%aaa6 = tic;
ii = (numx+indL+indE+indS)+indE+indS+1:(numx+indL+indE+indS)+indE+2*indS; %target cells

    
    %dVdt(:,i) = lambda-dt*V(i,:)-beta*V(i,:).*V(i-indS,:);
    %dVdt(:,i) = lambda-dt*V(i,:)-cellInfections-collisionsCC_t;
    
    %q = 1./(1+(IC50./(V(ii-indS))));
    %q = q(:);

    dVdt(ii,:) = lambda_S-dT_S.*V(ii,:)-(1-q).*cellInfections;

%fprintf('%30s %.5f %s\n', 'Timing', toc(aaa6), 'seconds');
%aaa7 = tic;
ii = (numx+indL+indE+indS)+indE+2*indS+1:(numx+indL+indE+indS)+indE+3*indS; %latent cells
    
    %dVdt(:,i) = eta*beta*V(i-indS,:).*V(i-indS-indS,:)-dl*V(i,:)-a*V(i,:);
    %MultConcCollisions = V(ii-indS,:).*V(ii-indS-indE-indS-indE-indG-indS-indS,:);
    %collisionsVC_t= 4*pi*(Rc+Rv)*(D_S)*MultConcCollisions;
    %cellInfections = perCollisionsInfection*collisionsVC_t;
    %q = 1./(1+(IC50./(V(ii-2*indS))));
    %q = q(:);

    dVdt(ii,:) = lambda_L*V(ii,:)+eta.*(1-q).*cellInfections-dL_S.*V(ii,:)-a.*V(ii,:);
%fprintf('%30s %.5f %s\n', 'Timing', toc(aaa7), 'seconds');
%aaa8 = tic;
ii = (numx+indL+indE+indS)+indE+3*indS+1:(numx+indL+indE+indS)+indE+4*indS; %infected cells
    %dVdt(i) = beta*V(i-ind2-ind2).*V(i-ind2-ind2-ind2)-del*V(i);
    %if t<tau_0
        %dVdt(:, i) = (1-eta)*beta*V(i-indS-indS,:).*V(i-indS-indS-indS,:)-del*V(i,:)+a*V(i-indS,:);
        %dVdt(:, i) = (1-eta)*cellInfections+collisionsCC_t-del*V(i,:)+a*V(i-1,:);
    %MultConcCollisions = V(ii-indS-indS,:).*V(ii-indS-indE-indS-indE-indG-indS-indS-indS,:);
    %collisionsVC_t= 4*pi*(Rc+Rv)*(D_S)*MultConcCollisions;
    %cellInfections = perCollisionsInfection*collisionsVC_t;
    %q = 1./(1+(IC50./(V(ii-3*indS))));
    %q = q(:);

    dVdt(ii,:) = (1-eta).*(1-q).*cellInfections-dI_S.*V(ii,:)+a.*V(ii-indS,:);
    
%% blood compartment

%%% TFV in blood
ii = (numx+indL+indE+indS)+indE+4*indS+1;


    %volumeCell = 4/3*pi*(10^(-4)^2/2)^3;

    %tempT = linspace(0, t, 100);
    %Dd_G,Dd_E,Dd_S,phid_GE,phid_ES,kd_D,kd_B,kd_L,Vdb
    % numx+indG+indE+1
    % numx+indG+indE+indS
    %CS = trapz(V(numx+indG+indE+1:numx+indG+indE+indS, :))/(indS-1);
    CS = trapz(V(numx+indL+indE+1:numx+indL+indE+indS, :))/(indS-1);
    %CS = mean(V(numx+indG+indE+1:numx+indG+indE+indS, :));
    Msbd = 2 * CS * W * L * h_S * kd_B;
    
    dVdt(ii,:) = (1/Vdb) * (Msbd - kd_L * V(ii,:) * Vdb);
    %kd_L * V(ii,:)
    %V(ii,:)

%%% TFV-DP in blood
ii = (numx+indL+indE+indS)+indE+4*indS+2;

    cellFractionBlood = cellV*V(ii+1,:)/Vvb;
    bracks = (V(ii-1,:).*cellFractionBlood-V(ii,:)/r);
    dVdt(ii,:) = Kon*bracks.*(bracks>0)-Koff*V(ii,:);
    
%%%%TENOFOVIR CORRECTION DUE TO TENOFOVIR DIPHOSPHATE CONVERSION IN
%%%%BLOOD

    dVdt((numx+indL+indE+indS)+indE+4*indS+1,:) = dVdt((numx+indL+indE+indS)+indE+4*indS+1,:) - dVdt((numx+indL+indE+indS)+indE+4*indS+2,:);

    qb = calcQ(V((numx+indL+indE+indS)+indE+4*indS+1,:));

%%% T in blood
ii = (numx+indL+indE+indS)+indE+4*indS+3;
    
    dVdt(ii,:) = lambda_B-(1-qb).*betaInf.*V(ii+3,:).*V(ii,:)-dT_B*V(ii,:);

%%% L in blood

ii = (numx+indL+indE+indS)+indE+4*indS+4;

    dVdt(ii,:) = lambda_L*V(ii,:)+eta.*(1-qb).*betaInf.*V(ii+2,:).*V(ii-1,:)-dL_B*V(ii,:)-a*V(ii,:);

%%% I in blood

ii = (numx+indL+indE+indS)+indE+4*indS+5;

    dVdt(ii,:) = (1-eta)*(1-qb).*betaInf.*V(ii+1).*V(ii-2,:)-dI_B*V(ii,:)+a*V(ii-1,:);
%     dIdt = (1-eta)*beta*BB(4)*BB(1)-d_I*BB(3)+a*BB(2);
%     dVdt = (1/(Vvb))*(Msb_V-kv_L*BB(4)+rho*BB(3));

%%% V in blood
ii = (numx+indL+indE+indS)+indE+4*indS+6;

    %tempT = linspace(0, t, 100);
    %Dd_G,Dd_E,Dd_S,phid_GE,phid_ES,kd_D,kd_B,kd_L,Vdb
    % numx+indG+indE+1
    % numx+indG+indE+indS
    %CS = trapz(V(numx+indG+indE+1:numx+indG+indE+indS, :))/(indS-1);
    VS = trapz(V((indL+indE+2):numx-1,:))/(indS-1);
    %CS = mean(V(numx+indG+indE+1:numx+indG+indE+indS, :));
    Msb = 2 * VS * W * L * h_S * kv_B;
    
    dVdt(ii,:) = (1/Vvb) * (Msb - kv_L * V(ii,:) * Vvb);

DVDT = dVdt;

end

% function q = PotCalc(C,f,maxPot,minPot,maxConc,minConc)
% 
%  
%     %EVALUATING CONCENTRATION TFV-DP IN TISSUE USING FIT 
%     evalPot = f.a0 + f.a1*cos(C*f.w) + f.b1*sin(C*f.w) + ...
%           f.a2*cos(2*C*f.w) + f.b2*sin(2*C*f.w) + ...
%           f.a3*cos(3*C*f.w) + f.b3*sin(3*C*f.w) + ...
%           f.a4*cos(4*C*f.w) + f.b4*sin(4*C*f.w) + ...
%           f.a5*cos(5*C*f.w) + f.b5*sin(5*C*f.w) + ...
%           f.a6*cos(6*C*f.w) + f.b6*sin(6*C*f.w) + ...
%           f.a7*cos(7*C*f.w) + f.b7*sin(7*C*f.w) + ...
%           f.a8*cos(8*C*f.w) + f.b8*sin(8*C*f.w);
%     
%     pot = minPot*(C<=minConc)+maxPot*(C>=maxConc)+evalPot.*(C>minConc & C<maxConc);
%     pot(C<=0) = minPot;
% 
%     q = 1-((pot-maxPot)/(minPot-maxPot));
%     q = q(:);
% 
% end


function S = calcsparsity(numx,indL,indE,indS)
totalSize = (2*numx)+(indS*4)+indE+6;
B = ones(totalSize,totalSize);
B(numx+1:2*numx+indE+4*indS,1:indL+indE) = 0; %a
B(1:indL+indE,numx+1:2*numx+indE+4*indS) = 0; %a
B(numx+1:numx+indL,2*numx+1:2*numx+indE+4*indS) = 0; %b
B(2*numx+1:2*numx+indE+4*indS,numx+1:numx+indL) = 0; %b
B(2*numx+indE+1:2*numx+indE+4*indS,numx+indL+1:numx+indL+indE) = 0; %c
B(numx+indL+1:numx+indL+indE,2*numx+indE+1:2*numx+indE+4*indS) = 0; %c
B(2*numx+indE+indS+1:2*numx+indE+4*indS,numx+indL+1:numx+indL+indE)=0; %d
B(numx+indL+1:numx+indL+indE,2*numx+indE+indS+1:2*numx+indE+4*indS)=0; %d
B(2*numx+indE+1:2*numx+indE+4*indS,2*numx+1:2*numx+indE)=0; %e
B(2*numx+1:2*numx+indE,2*numx+indE+1:2*numx+indE+4*indS)=0; %e
B(2*numx+1:2*numx+indE,numx+indL+indE+1:2*numx)=0; %f
B(numx+indL+indE+1:2*numx,2*numx+1:2*numx+indE)=0; %f

%B(numx+1:2*numx+indE+indS,indG+indE+1:numx) = 0;
%B(indG+indE+1:numx,numx+1:2*numx+indE+indS) = 0;

S = B;

end
