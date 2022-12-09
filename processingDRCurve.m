function [f,maxPot,minPot,maxConc,minConc] = processingDRCurve(T_0)

%CALCULATING MICROBICIDE EFFICACY Q

    %Original data from 
    TFVDPPotency = load('TFVDP_Potency.txt');
    smoothTFVDPPotency = smooth(TFVDPPotency(:,1),TFVDPPotency(:,2),200,'rloess');
    % ConcX = log10(10.^(TFVDPPotency(:,1))*(10^(6)/180)); %4.9917*10^3 10^6cells / mL tissue  = 4.9917 10^6cells / mg tissue. 
    %                                                %Assuming: density tissue = 1000mg/mL. Conversion from Schwartz et al.
    %                                                %fmol/10^6 cells to fmol/mg tissue
    ConcX = log10(10.^(TFVDPPotency(:,1))*(T_0/(10^(6))));
    %ConcX = log10(10.^(TFVDPPotency(:,1))*4.9917*10^(-2));
    
    %MAKING INITIAL FIT

    f = fit(ConcX,smoothTFVDPPotency,'fourier8');
    fvals = f.a0 + f.a1*cos(ConcX*f.w) + f.b1*sin(ConcX*f.w) + ... %Fit for concentration in units of log10 fmol/mg tissue
            f.a2*cos(2*ConcX*f.w) + f.b2*sin(2*ConcX*f.w) + ...
            f.a3*cos(3*ConcX*f.w) + f.b3*sin(3*ConcX*f.w) + ...
            f.a4*cos(4*ConcX*f.w) + f.b4*sin(4*ConcX*f.w) + ...
            f.a5*cos(5*ConcX*f.w) + f.b5*sin(5*ConcX*f.w) + ...
            f.a6*cos(6*ConcX*f.w) + f.b6*sin(6*ConcX*f.w) + ...
            f.a7*cos(7*ConcX*f.w) + f.b7*sin(7*ConcX*f.w) + ...
            f.a8*cos(8*ConcX*f.w) + f.b8*sin(8*ConcX*f.w);
    maxPot = min(fvals);
    minPot = max(fvals);
    maxConc = max(ConcX);
    minConc = min(ConcX);
    
end