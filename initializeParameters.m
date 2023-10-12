
function [V, T, I, U, D, Geometry] = initializeParameters(params)
    V = initializeVirus(params);
    T = initializeTargetCells(params);
    I = initializeInfectedCells(params);
    U = initializeChemokine(params);
    D = initializeDrug(params);
    Geometry = initializeGeometry(params);
end

function D = initializeDrug(params)
    D = Drug;
    D.initialLumenConcentration = params.C_G0;
    D.diffCoeffLumen = 6*10^(-6);
    D.diffCoeffEpithelium = 7*10^(-8);
    D.diffCoeffStroma = 4*10^(-7);
    D.partitionLumenEpthelium = 0.75;
    D.partitionEpitheliumStroma = 1;
    D.lossDilutionInLumen = 1.22/3600;
    D.lossStromaToBlood = 0.119/3600;
    D.lossClearanceInBlood = 1.41/3600;
    D.volumeOfDistribution = 75*1000;
    D.rateActivation = log(2)/3600;
    D.rateDeactivation = log(2)/(7*24*3600);
    D.ratioActivationFromDeactivated = 0.1;
    D = D.setSymbols();
end

function V = initializeVirus(params)
    V = Virus;
    V.diffCoeffLumen = 1.27*10^(-8); %D_G
    V.diffCoeffTissue = 1*10^(-9); %D_T
    V.diffCoeffStroma = V.diffCoeffTissue; %D_S
    V.diffCoeffEpithelium = 1*10^(-10); %D_E
    V.lossDilutionInLumen = 0.3/3600; %k_D
    V.lossTissueToBlood = params.k_B; %k_B %%% MAY NEED TO FIND THIS
    V.lossClearanceInBlood = params.k_L; %k_L
    V.partitionLumenTissue = 0.3; %phi_GE %Hope et al, 2015. In Macaques, and Humans
    V.volumeOfDistribution = 5*1000; %volume of blood in ml
    %V.volumeOfDistribution = 1e15*1000; %same as TFV jings paper
    V.initialLumenConcentration = params.V_0; %V_0 in virions/ml
    V.radius = 100*10^(-7)/2; %in cm
    V.ratioCollisionsCausingInfection = 1/10; % 10%
    %V.infectivity = (6.5*10^(-7))/(24*3600); %perelson paper
    %V.infectivity = (1e-5)/(24*3600); % perelson annual review paper 2016
    V.infectivity = params.beta;
    V = V.setSymbols();
end

function T = initializeTargetCells(params)
    T = TCell;
    %T.initialTissueConcentration = 10^(4)*1.088*10^(3)*1/100; %T_0 *1.088*10^(3)*1/100 cells/ml
    %T.initialTissueConcentration = 5e5;
    T.initialTissueConcentration = params.T0_T;
    %T.initialBloodConcentration = 2e6; %T_0 PBMC
    %T.initialBloodConcentration = 1e6; %T_0 PBMC Between 500-1500cells/mm3 AIDS CD4+ Count, Guzman, 2022
    T.initialBloodConcentration = params.T0_B;
    T.volumeFractionTissue = 0.1; %phiDP_S = 0.1; jing 2D 
    T.deathRateTissue = 0.01/(24*3600); %s^(-1)
    T.deathRateBlood = 0.01/(24*3600); %s^(-1)
    %T.deathRateTissue = 0.1/(24*3600); %s^(-1) % perelson annual review paper 2016
    %T.deathRateBlood = 0.1/(24*3600); %s^(-1) % perelson annual review paper 2016
    T.productionRateTissue = T.deathRateTissue * T.initialTissueConcentration;
    T.productionRateBlood = T.deathRateBlood * T.initialBloodConcentration;
    T.radius = 10*10^(-4)/2; %diameter = 10um
    T.diffCoeffTissue = (3*10^(-4))/(3600*(10)^2); %3e-4mm2/hr haugh et al, 2006
    T.uptakeRateFromBlood = 1/3600; %this is in units of (concentration of chemokine)^-1 (s)^-1 THIS IS UNKNOWN
    T.maxChemotaxSpeed = 0.1/(10*3600); %0.1 mm/hr haugh et al 2006
    T.sensitivityToChemokine = 5; %HAVE NO IDEA WHAT TO PUT HERE
    T = T.setSymbols();
end

function I = initializeInfectedCells(params)
    I = ICell;
    I.productionRateVirus = params.rho; %rho
    I.deathRateTissue = 0.39/(24*3600); %dI stroma
    I.deathRateBlood = 0.39/(24*3600);
    %I.deathRateTissue = 0.5/(24*3600); % perelson annual review paper 2016 
    %I.deathRateBlood = 0.5/(25*3600);
    I.diffCoeffTissue = (3*10^(-4))/(3600*(10)^2); %3e-4 haugh et al, 2006
    I.lossTissueToBlood = 0.1/3600; %1/hr approximated
    %I.cellToCellInfection = (1e-3)/(24*3600); %s^(-1) per infected cell. Viral transmission
    %I.cellToCellInfection = (2e-4)/(24*3600); % USING THIS AS DEFAULT
    I.cellToCellInfection = params.w;
    I = I.setSymbols();
end

function U = initializeChemokine(params)
    U = Chemokine;
    U.lossRateBinding = 1/3600; %1/hr haugh et al
    U.lossRateDegradation = 0.1/3600; %0.1/hr deterministic model of dermal wound invasion, haugh
    U.secretionRate = 1/3600; %1/3600 - 1/hr haugh et al.
    U.diffCoeffTissue = 0.01/(3600*(10)^2); %0.01 mm2/hr = 2.78e-8 cm2/s, haugh et al 2006
    U = U.setSymbols();
end

function Geometry = initializeGeometry(params)
    Geometry = Vagina;
    Geometry.thicknessLumen = 0.04; %h_G cm
    Geometry.thicknessEpithelium = params.h_E; %h_E cm
    Geometry.thicknessStroma = 0.28; %h_S cm
    Geometry.thicknessTissue = Geometry.thicknessEpithelium + Geometry.thicknessStroma;
    Geometry.width = 3.35; % w cm
    Geometry.length = 13; %l cm
    Geometry.volumeFractionCellsEpithelium = 0.95; %phi_e
    Geometry.volumeFractionCellsStroma = 0.1; %phi_s
    Geometry = Geometry.setSymbols();
end

