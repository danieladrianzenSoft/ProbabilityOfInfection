function newparams = samplingNormalPOI_BinaryCR(filename,N,paramMeans,overwrite)

%% Viral Load (V_0)

if isKey(paramMeans,"V_0") == 1
    mu = paramMeans("V_0");
else
    mu = 1*10^(4);
end
sigma = (1*10^(4)/(1*10^(4)))*mu;

lognormalMu = log(mu^2/(sqrt(mu^2 + sigma^2)));
lognormalSigma = log(1+(sigma^2/mu^2));

pd = makedist('lognormal',lognormalMu,lognormalSigma);
t=truncate(pd,10,inf);
V_0 = random(t,N,1);

%% Epithelium Thickness (h_E)

% input: cm

if isKey(paramMeans,"h_E") == 1
    mu = paramMeans("h_E");
else
    mu = 0.020; %cm
end
sigma = (0.002*sqrt(15)/0.02)*mu;

lognormalMu = log(mu^2/(sqrt(mu^2 + sigma^2)));
lognormalSigma = log(1+(sigma^2/mu^2));

pd = makedist('lognormal',lognormalMu,lognormalSigma);
t=truncate(pd,0.001,inf);
h_E = random(t,N,1);
%influence of hormonal contraceptives on the immune cells and thickness of
%vaginal epithelium, ildgruben, 2003.
%epithelial cell layer thickness and immune cell populations in the normal
%human vagina at different stages of the menstrual cycle, dorothy patton,
%1999.

%% Production rate of V from I (rho = pi)

% input: virs/day
% convert to: virs/s

if isKey(paramMeans,"rho") == 1
    mu = paramMeans("rho")/(24*3600);
else
    mu = 850/(24*3600);
end
%mu = (5*10^3/(24*3600)); %virion / s
%sigma = (1000/(24*3600));
sigma = (1*10^3/850) * mu;

lognormalMu = log(mu^2/(sqrt(mu^2 + sigma^2)));
lognormalSigma = log(1+(sigma^2/mu^2));

pd = makedist('lognormal',lognormalMu,lognormalSigma);
t=truncate(pd,0.1/(24*3600),inf);
rho = random(t,N,1);

%modeling plasma virus concentration during primary infection, stafford,
%2000
%Current estimates for HIV production imply rapid viral clearance in
%lymphoid tissues. De boer, 2010.

%% clearance rate from tissue to blood(k_B)
% input: /day
% convert to: /s

if isKey(paramMeans,"k_B") == 1
    mu = paramMeans("k_B")/(24*3600);
else
    mu = 6/(24*3600);
end
sigma = (6 / 6) * mu;

lognormalMu = log(mu^2/(sqrt(mu^2 + sigma^2)));
lognormalSigma = log(1+(sigma^2/mu^2));

pd = makedist('lognormal',lognormalMu,lognormalSigma);
t=truncate(pd,0.05/(24*3600), inf);
k_B = random(t,N,1);

%% clearance rate from blood (k_L = c)

% input: /day
% convert to: /s
if isKey(paramMeans,"k_L") == 1
    mu = paramMeans("k_L")/(24*3600);
elseif isKey(paramMeans,"c") == 1
    mu = paramMeans("c")/(24*3600);
else
    mu = 20/(24*3600);
end
%mu = 20/(24*3600); 
sigma = (10 / 20) * mu;

lognormalMu = log(mu^2/(sqrt(mu^2 + sigma^2)));
lognormalSigma = log(1+(sigma^2/mu^2));

pd = makedist('lognormal',lognormalMu,lognormalSigma);
t=truncate(pd,3/(24*3600), inf);
k_L = random(t,N,1);

%Current estimates for HIV production imply rapid viral clearance in
%lymphoid tissues. De boer, 2010.


%% Infectivity (beta = k)

if isKey(paramMeans,"beta") == 1
    mu = paramMeans("beta")/(24*3600);
else
    mu = (0.65*10^(-6))/(24*3600); %cm^3 virion / s
end
%sigma = (1.4*10^(-6)/(24*3600));
sigma = (1.4*10^(-6) / (0.65*10^(-6))) * mu;

lognormalMu = log(mu^2/(sqrt(mu^2 + sigma^2)));
lognormalSigma = log(1+(sigma^2/mu^2));

pd = makedist('lognormal',lognormalMu,lognormalSigma);
t=truncate(pd,1*10^(-10)/(24*3600),inf);
beta = random(t,N,1);

%% Cell-to-cell infectivity (w)

if isKey(paramMeans,"w") == 1
    w = (paramMeans("w")/(24*3600))*ones(N,1);
else
    w = ((2e-4)/(24*3600))*ones(N,1);
end

%% T0 in Blood

if isKey(paramMeans,"T0_B") == 1
    mu = paramMeans("T0_B");
else
    mu = 1e4; %cells / cm^3
end
%sigma = (1.4*10^(-6)/(24*3600));
sigma = (3e3 / 1e4) * mu;

lognormalMu = log(mu^2/(sqrt(mu^2 + sigma^2)));
lognormalSigma = log(1+(sigma^2/mu^2));

pd = makedist('lognormal',lognormalMu,lognormalSigma);
t=truncate(pd,1e2,inf);
T0_B = random(t,N,1);

%% T0 in Tissue (Stroma)

if isKey(paramMeans,"T0_T") == 1
    mu = paramMeans("T0_T");
else
    mu = 1e3; %cells / cm^3
end
%sigma = (1.4*10^(-6)/(24*3600));
sigma = (3e2 / 1e3) * mu;

lognormalMu = log(mu^2/(sqrt(mu^2 + sigma^2)));
lognormalSigma = log(1+(sigma^2/mu^2));

pd = makedist('lognormal',lognormalMu,lognormalSigma);
t=truncate(pd,10,inf);
T0_T = random(t,N,1);
% %HIV Diffusion coefficient D_vS 
% 
% mu = 5*10^(-9); 
% sigma = 3*10^(-9);
% 
% pd = makedist('Normal',mu,sigma);
% t=truncate(pd,1*10^(-10),inf);
% D_vS = random(t,N,1);

% Surface Area
% input: cm2

if isKey(paramMeans,"SA") == 1
    SA = paramMeans("SA")*ones(N,1);
else
    SA = zeros(N,1);
end

% Volume of vaginal fluid
% input: cm3

if isKey(paramMeans,"Vf") == 1
    Vf = paramMeans("Vf")*ones(N,1);
else
    Vf = zeros(N,1);
end

% Volume of semen
% input: cm3

if isKey(paramMeans,"Vs") == 1
    Vs = paramMeans("Vs")*ones(N,1);
else
    Vs = zeros(N,1);
end

% Volume of drug gel 
% input: cm3

if isKey(paramMeans,"Vs") == 1
    Vg = paramMeans("Vg")*ones(N,1);
else
    Vg = zeros(N,1);
end

% Initial Drug Concentration

if isKey(paramMeans,"C0") == 1
    C0 = paramMeans("C0")*ones(N,1);
else
    C0 = zeros(N,1);
end

% Drug application vs HIV exposure delay T_VD 
% input: hrs, 
% convert to: s

if isKey(paramMeans,"T_VD") == 1
    T_VD = (paramMeans("T_VD")*3600)*ones(N,1); 
else
    T_VD = zeros(N,1);
end

% Menstrual cycle phase (0 = midcycle, 1 = follicular, 2 = luteal)

if isKey(paramMeans,"phase") == 1
    phase = paramMeans("phase")*ones(N,1);
else
    phase = zeros(N,1);
end

% Drug type (0 = none, 1 = gel, 2 = IVR)

if isKey(paramMeans,"drugType") == 1
    if (isKey(paramMeans,"C0") == 1 && paramMeans("C0") == 0 && paramMeans("drugType") ~= 0)
        drugType = 0*ones(N,1);
        fprintf('WARNING: Incongruence, Inputted C0 = %d, but drugType = %d. Using drugType = 0\n', paramMeans("C0"), paramMeans("drugType"));
    else
        drugType = paramMeans("drugType")*ones(N,1);
    end
else
    drugType = zeros(N,1);
end


% T_VD = paramMeans("T_VD")*ones(N,1);

% if isKey(paramMeans,"isRandomT_VD") == 1 && isKey(paramMeans,"T_VD") == 1
%     if paramMeans("isRandomT_VD") == 0
%         T_VD = paramMeans("T_VD")*ones(N,1);
%     else
%         % sample from uniform distribution between -6hrs and 6hrs
%         T_VD = -6 + (6+6)*rand(N,1);
%     end
% else
%     T_VD = zeros(N,1);
% end


newparams = table(V_0,SA,h_E,rho,k_B,k_L,beta,w,T0_B,T0_T,Vf,Vs,Vg,C0,T_VD,phase,drugType);

if overwrite == 1
    writetable(newparams,filename,'WriteMode','overwritesheet',...
    'WriteRowNames',true)
else
    writetable(newparams,filename,'WriteMode','Append',...
    'WriteVariableNames',false,'WriteRowNames',true) 
end

