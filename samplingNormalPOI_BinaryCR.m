function newparams = samplingNormalPOI_BinaryCR(filename,N,C0,DVdelay,overwrite)

%% Viral Load (V_0)

mu = 1*10^(4);
sigma = 1*10^(4);

lognormalMu = log(mu^2/(sqrt(mu^2 + sigma^2)));
lognormalSigma = log(1+(sigma^2/mu^2));

pd = makedist('lognormal',lognormalMu,lognormalSigma);
t=truncate(pd,400,inf);
V_0 = random(t,N,1);

%% Epithelium Thickness (h_E)

mu = 0.020; %cm
sigma = 0.002*sqrt(15);

lognormalMu = log(mu^2/(sqrt(mu^2 + sigma^2)));
lognormalSigma = log(1+(sigma^2/mu^2));

pd = makedist('lognormal',lognormalMu,lognormalSigma);
t=truncate(pd,0.005,inf);
h_E = random(t,N,1);
%influence of hormonal contraceptives on the immune cells and thickness of
%vaginal epithelium, ildgruben, 2003.
%epithelial cell layer thickness and immune cell populations in the normal
%human vagina at different stages of the menstrual cycle, dorothy patton,
%1999.

%% Production rate of V from I (rho = pi)

%mu = (5*10^3/(24*3600)); %virion / s
mu = 850/(24*3600);
%sigma = (1000/(24*3600));
sigma = (1*10^3/(24*3600));

lognormalMu = log(mu^2/(sqrt(mu^2 + sigma^2)));
lognormalSigma = log(1+(sigma^2/mu^2));

pd = makedist('lognormal',lognormalMu,lognormalSigma);
t=truncate(pd,5/(24*3600),inf);
rho = random(t,N,1);

%modeling plasma virus concentration during primary infection, stafford,
%2000
%Current estimates for HIV production imply rapid viral clearance in
%lymphoid tissues. De boer, 2010.

%% clearance rate from tissue to blood(k_B = c)

mu = 6/(24*3600); 
sigma = 2/(24*3600);

lognormalMu = log(mu^2/(sqrt(mu^2 + sigma^2)));
lognormalSigma = log(1+(sigma^2/mu^2));

pd = makedist('lognormal',lognormalMu,lognormalSigma);
t=truncate(pd,1/(24*3600), inf);
k_B = random(t,N,1);

%% clearance rate from blood (k_L = c)

mu = 20/(24*3600); 
sigma = 8/(24*3600);

lognormalMu = log(mu^2/(sqrt(mu^2 + sigma^2)));
lognormalSigma = log(1+(sigma^2/mu^2));

pd = makedist('lognormal',lognormalMu,lognormalSigma);
t=truncate(pd,3/(24*3600), inf);
k_L = random(t,N,1);

%Current estimates for HIV production imply rapid viral clearance in
%lymphoid tissues. De boer, 2010.


%% Infectivity (beta = k)

mu = ((0.65*10^(-6))/(24*3600))/10; %cm^3 virion / s
sigma = (1.4*10^(-6)/(24*3600))/5;

pd = makedist('Normal',mu,sigma);
t=truncate(pd,1*10^(-9)/(24*3600),5*10^(-6)/(24*3600));
beta = random(t,N,1);

% %HIV Diffusion coefficient D_vS 
% 
% mu = 5*10^(-9); 
% sigma = 3*10^(-9);
% 
% pd = makedist('Normal',mu,sigma);
% t=truncate(pd,1*10^(-10),inf);
% D_vS = random(t,N,1);

%Initial Drug Concentration

C_G0 = C0*ones(N,1);

%Drug application vs HIV exposure delay T_VD 

T_VD = DVdelay*ones(N,1);


newparams = table(V_0,h_E,rho,k_B,k_L,beta,C_G0,T_VD);

if overwrite == 1
    writetable(newparams,filename,'WriteMode','overwritesheet',...
    'WriteRowNames',true)
else
    writetable(newparams,filename,'WriteMode','Append',...
    'WriteVariableNames',false,'WriteRowNames',true) 
end

