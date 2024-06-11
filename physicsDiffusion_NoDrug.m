function DVDT = physicsDiffusion_NoDrug(t,sol,IC,Geometry,Mesh,V,T,I,D,U,includeOutputs)
    %Here we add in cell movement via diffusion only, chemokine is
    %still not included. Collision theory prevails in the tissue
    %infectivity in the blood. 
    %4 compartments: lumen, epithelium, stroma, blood.

    %NOTE: TO UNDERSTAND THESE EQUATIONS, LOOK UP FOR DIFFERENCE EQUATIONS: CENTRAL, BACKWARD
    %AND FORWARD.
    
    dsoldt = zeros(length(IC),size(sol,2));
    dsoldx = zeros(length(IC),size(sol,2));
    d2soldx2 = zeros(length(IC),size(sol,2));
    
    dxL = Geometry.h_L / (Mesh.numL - 1);
    dxE = Geometry.h_E / (Mesh.numE - 1);
    dxS = Geometry.h_S / (Mesh.numS - 1);

    %INDICES
    V_LIndices = 1 : Mesh.numL;
    V_intf_LE = Mesh.numL;
    V_EIndices = Mesh.numL + 1 : Mesh.numL + Mesh.numE;
    V_intf_ES = Mesh.numL + Mesh.numE;
    V_SIndices = Mesh.numL + Mesh.numE + 1 : Mesh.numX;
    T_SIndices =  Mesh.numX + 1 : Mesh.numX + Mesh.numS;
    I_SIndices =  Mesh.numX + Mesh.numS + 1 : Mesh.numX + 2 * Mesh.numS;
    V_BIndices =  Mesh.numX + 2 * Mesh.numS + 1;
    T_BIndices = Mesh.numX + 2 * Mesh.numS + 2;
    I_BIndices = Mesh.numX + 2 * Mesh.numS + 3;

    if (includeOutputs == 1)
        VirusCellStroma = Mesh.numX + 2 * Mesh.numS + 4;
        CellCellStroma = Mesh.numX + 2 * Mesh.numS + 5;
        VirusCellBlood = Mesh.numX + 2 * Mesh.numS + 6;
        CellCellBlood = Mesh.numX + 2 * Mesh.numS + 7;
    end

    %HIV INTERFACE
    V_intf1a = (V.DV_L.*sol(V_intf_LE,:)+(V.DV_E.*sol(V_intf_LE+1,:)))./((V.DV_E.*V.phi_LT)+V.DV_L); %V at lumen/epithelium interface
    V_intf1b = (V.DV_L.*sol(V_intf_LE,:)+(V.DV_E.*sol(V_intf_LE+1,:)))./((V.DV_L./V.phi_LT)+V.DV_E);
    V_intf2a = (V.DV_E.*sol(V_intf_ES,:)+(V.DV_S.*sol(V_intf_ES+1,:)))./(V.DV_S+V.DV_E); % V at epithelium/tissue interface, assuming part.coeff of 1
    V_intf2b = (V.DV_E.*sol(V_intf_ES,:)+(V.DV_S.*sol(V_intf_ES+1,:)))./(V.DV_E+V.DV_S);
    

    %% %%%%%%%%% TISSUE %%%%%%%%%% %%
    
    %%%%% VIRUS %%%%%%
    
    ii = V_LIndices(1); %y direction BC zero flux (dcdx=0)
        d2soldx2(ii,:) = 2 * (sol(ii+1,:) - sol(ii,:))./(dxL.^2);
        dsoldt(ii,:) = V.DV_L .* d2soldx2(ii,:) - (V.kD.*sol(ii,:));
        %dVdt(i) = 2*D_G.*(V(i+1)-V(i))./(dx.^2);
        
    ii = V_LIndices(2:end-1); %in lumen
        %dVdx(ii,:) = (V(ii+1,:)-V(ii,:))./dxG;
        d2soldx2(ii,:) = (sol(ii+1,:) - 2.*sol(ii,:) + sol(ii-1,:))./(dxL.^2);
        dsoldt(ii,:) = V.DV_L.*d2soldx2(ii,:) - (V.kD.*sol(ii,:));
    
    ii = V_LIndices(end); %right before interface - lumen
        %dVdx(ii,:) = (V_intf1a - V(ii,:))./dxG;
        d2soldx2(ii,:) = (V_intf1a - 2.*sol(ii,:) + sol(ii-1,:))./(dxL.^2);
        dsoldt(ii,:) = V.DV_L.*d2soldx2(ii,:) - (V.kD.*sol(ii,:));
    
    ii = V_EIndices(1); %right after interface - in epithelium
        %dVdx(ii,:) = (V(ii+1,:)-V(ii,:))./dxE;
        d2soldx2(ii,:) = (sol(ii+1,:) - 2.*sol(ii,:) + (V_intf1b))./(dxE.^2);
        %dsoldt(ii,:) = V.DV_E.*d2soldx2(ii,:) - (V.kL).*sol(ii,:);
        dsoldt(ii,:) = V.DV_E.*d2soldx2(ii,:);
    
    ii = V_EIndices(2:end-1); %in epithelium
        %dVdx(ii,:) = (V(ii+1,:)-V(ii,:))./dxE;
        d2soldx2(ii,:) = (sol(ii+1,:) - 2.*sol(ii,:) + sol(ii-1,:))./(dxE.^2);
        %dsoldt(ii,:) = V.DV_E.*d2soldx2(ii,:) - (V.kL).*sol(ii,:);
        dsoldt(ii,:) = V.DV_E.*d2soldx2(ii,:);
    
    ii = V_EIndices(end); %right before interface - in epithelium
        %dVdx(ii,:) = (V_intf2a-V(ii,:))./dxE;
        d2soldx2(ii,:) = (V_intf2a - 2.*sol(ii,:) + sol(ii-1,:))./(dxE.^2);
        %dsoldt(ii,:) = V.DV_E*d2soldx2(ii,:) - (V.kL).*sol(ii,:); %BC zero flux
        dsoldt(ii,:) = V.DV_E.*d2soldx2(ii,:);

    ii = V_SIndices(1); %right after interface - in stroma
        %dVdx(ii,:) = (V(ii+1,:)-V(ii,:))./dxE;
        d2soldx2(ii,:) = (sol(ii+1,:) - 2.*sol(ii,:) + (V_intf2b))./(dxS.^2);
        dsoldt(ii,:) = V.DV_S.*d2soldx2(ii,:) - (V.kB).*sol(ii,:) + I.rho*sol(I_SIndices(1),:);
    
    ii = V_SIndices(2:end-1); %in stroma
        %dVdx(ii,:) = (V(ii+1,:)-V(ii,:))./dxE;
        d2soldx2(ii,:) = (sol(ii+1,:) - 2.*sol(ii,:) + sol(ii-1,:))./(dxS.^2);
        dsoldt(ii,:) = V.DV_S.*d2soldx2(ii,:) - (V.kB).*sol(ii,:) + I.rho*sol(I_SIndices(2:end-1),:);
    
    ii = V_SIndices(end); %end of stroma
        %dVdx(ii,:) = (V_intf2a-V(ii,:))./dxE;
        d2soldx2(ii,:) = 2*(sol(ii-1,:)-sol(ii,:))./(dxS.^2);
        dsoldt(ii,:) = V.DV_S*d2soldx2(ii,:) - (V.kB).*sol(ii,:) + I.rho*sol(I_SIndices(end),:); %BC zero flux            

    %%%%% CELLS %%%%%%
    
    %zeta_T = 4 * pi * (T.r + V.r) .* (T.DT_T + V.DV_T); 
    % No cell diffusion
    %zeta_T = 4 * pi * (T.r + V.r) .* (V.DV_S);

    % Uavg = trapz(sol(U_TIndices, :)) / (Mesh.numT - 1);
    
    
    ii = T_SIndices(1); %target cells
        
        %Vb_tissue = 2 * Geometry.W * Geometry.L * Geometry.h_T;
        %dsoldx(ii,:) = (sol(ii+1,:)-sol(ii-1,:))./(2*dxT);
        d2soldx2(ii,:) = 2.*(sol(ii+1,:)-sol(ii,:))./(dxS.^2);
        dsoldt(ii,:) = ...
                     T.DT_T.*d2soldx2(ii,:) ... %diffusion
                     + (T.lambda_T - T.dT_T .* sol(ii,:)) ... % proliferation rate
                     - I.w .* sol(I_SIndices(1),:) .* sol(ii,:)... % loss due to cell-cell infections
                     - V.beta .* sol(V_SIndices(1), :) .* sol(ii,:); % loss due to infections, collision theory
                     ...- V.s .* zeta_T .* sol(V_SIndices(1), :) .* sol(ii,:); % loss due to infections, collision theory
                     
                     %... T.Stax .* T.Sigma .* d2soldx2(U_TIndices(1),:) .* sol(ii,:); %chemotaxis with du/dx and dt/dx = 0 at x = hL
    
                     %Uavg .* sol(T_BIndices,:) .* ((Vb_tissue / V.Vb) .* T.Ku); %uptake from blood due to U(x,t)
                     %Uavg .* sol(T_BIndices,:) .* (T.Ku); %uptake from blood due to U(x,t)
                     %V.s .* zeta_T .* sol(Mesh.numL + 1, :) .* sol(ii,:); % loss due to infections
    
    ii = T_SIndices(2:end-1); %target cells
    
    %     zeta_T(2:end-1) = 4 * pi * (T.r + V.r) * (T.DT_T + V.DV_T) .* ...
    %     sol((Mesh.numL+2 : Mesh.numX-1),:) .* sol(ii,:); %collisions per volume per time
    
        dsoldx(ii,:) = (sol(ii+1,:)-sol(ii-1,:))./(2*dxS);
        d2soldx2(ii,:) = (sol(ii+1,:)-2.*sol(ii,:)+sol(ii-1,:))./(dxS.^2);
        dsoldt(ii,:) = ... 
                     T.DT_T.*d2soldx2(ii,:) ... %diffusion
                     + (T.lambda_T - T.dT_T .* sol(ii,:)) ... % proliferation rate,
                     - I.w .* sol(I_SIndices(2:end-1),:) .* sol(ii,:)... % loss due to cell-cell infections
                     - V.beta .* sol(V_SIndices(2:end-1), :) .* sol(ii,:);
                     ...- V.s .* zeta_T .* sol(V_SIndices(2:end-1),:).* sol(ii,:); % loss due to infections, collision theory


                    ...T.Stax .* T.Sigma .* (dsoldx(U_TIndices(2:end-1), :) .* sol(ii, :) + dsoldx(ii,:) .* sol(U_TIndices(2:end-1), :)); % chemotaxis first derivative u
                     
                    %Uavg .* sol(T_BIndices,:) .* ((Vb_tissue / V.Vb) .* T.Ku) + ... % uptake from blood due to U(x,t)
                     % T.Stax .* T.Sigma .* (d2soldx2(U_TIndices(2:end-1),:) .* sol(ii,:) + dsoldx(U_TIndices(2:end-1),:)) .* dsoldx(ii,:); %chemotaxis second derivative u
                     % V.s .* zeta_T .* sol((Mesh.numL+2):(Mesh.numL+Mesh.numT)-1,:) .* sol(ii,:); % loss due to infections
                     %V.s .* zeta_T .* sol((Mesh.numL+2):(Mesh.numL+Mesh.numT)-1,:) .* sol(ii,:); % loss due to infections
                     %Uavg .* sol(T_BIndices,:) .* T.Ku; %uptake from blood due to U(x,t)
        
    ii = T_SIndices(end); % target cells
    
        % dsoldx(ii,:) = (sol(ii+1,:)-sol(ii-1,:))./(2*dxT);
        d2soldx2(ii,:) = 2*(sol(ii-1,:)-sol(ii,:))./(dxS.^2);
        dsoldt(ii,:) = ...
                     T.DT_T.*d2soldx2(ii,:) ... %diffusion
                     + (T.lambda_T - T.dT_T .* sol(ii,:)) ... % proliferation rate
                     - I.w .* sol(I_SIndices(end),:) .* sol(ii,:) ... % loss due to cell-cell infections
                     - V.beta .* sol(V_SIndices(end), :) .* sol(ii,:);
                     ...- V.s .* zeta_T .* sol(V_SIndices(end),:).* sol(ii,:); % loss due to infections, collision theory
                     ... T.Stax .* T.Sigma .* d2soldx2(U_TIndices(end),:) .* sol(ii,:); %chemotaxis with du/dx and dt/dx = 0 at x = hL
    
                     %Uavg .* sol(T_BIndices,:) .* ((Vb_tissue / V.Vb) .* T.Ku); %uptake from blood due to U(x,t)
                     %V.s .* zeta_T .* sol(Mesh.numL + Mesh.numT) .* sol(ii,:); % loss due to infections 
                     %Uavg .* sol(T_BIndices,:) .* T.Ku; %updake from blood due to U(x,t)
    
    ii = I_SIndices(1); %infected cells
    
        d2soldx2(ii,:) = 2.*(sol(ii+1,:)-sol(ii,:))./(dxS.^2);
        dsoldt(ii,:) = ...
                     I.DI_T .* d2soldx2(ii,:) ... %diffusion
                     ...+ V.s.*zeta_T .* sol(V_SIndices(1), :) .* sol(T_SIndices(1),:) ... %increase due to infections, collision theory
                     + V.beta .* sol(V_SIndices(1), :) .* sol(T_SIndices(1),:) ...
                     + I.w .* sol(ii,:) .* sol(T_SIndices(1),:) ...
                     - I.dI_T * sol(ii,:); %death rate
                     % I.dI_T * sol(ii,:) - ... %death rate
                     % I.kI_TB * sol(ii,:); %loss to blood.
                
    ii = I_SIndices(2:end-1); %infected cells
    
        d2soldx2(ii,:) = (sol(ii+1,:)-2.*sol(ii,:)+sol(ii-1,:))./(dxS.^2);
        dsoldt(ii,:) = ... 
                     I.DI_T .* d2soldx2(ii,:) ... %diffusion
                     ...+ V.s.*zeta_T .* sol(V_SIndices(2:end-1),:) .* sol(T_SIndices(2:end-1),:) - ... %increase due to infections, collision theory
                     + V.beta .* sol(V_SIndices(2:end-1), :) .* sol(T_SIndices(2:end-1),:)...
                     + I.w .* sol(ii,:) .* sol(T_SIndices(2:end-1),:) ...
                     - I.dI_T * sol(ii,:); %death rate 
                     %I.dI_T * sol(ii,:) - ... %death rate
                     %I.kI_TB * sol(ii,:); %loss to blood.
    
        
    ii = I_SIndices(end); %infected cells
    
        d2soldx2(ii,:) = 2*(sol(ii-1,:)-sol(ii,:))./(dxS.^2);
        dsoldt(ii,:) = ...
                     I.DI_T.*d2soldx2(ii,:) ... %diffusion
                     ...+ V.s.*zeta_T .* sol(V_SIndices(end), :) .* sol(T_SIndices(end),:) ... %increase due to infections, collision theory
                     + V.beta .* sol(V_SIndices(end), :) .* sol(T_SIndices(end),:) ...
                     + I.w .* sol(ii,:) .* sol(T_SIndices(end),:) ...
                     - I.dI_T * sol(ii,:); %death rate
                     %I.dI_T * sol(ii,:) - ... %death rate
                     %I.kI_TB * sol(ii,:); %loss to blood.
    
    
    %% %%%%%%%%% BLOOD %%%%%%%%%% %%
    
    % zeta_B = 4 * pi * (T.r + V.r) * (T.DT_T + V.DV_T) .* ...
    %         sol(V_BIndices, :) .* sol(T_BIndices, :);
    
    %%% V in blood
    ii = V_BIndices;
    
        %VS = trapz(sol(Mesh.numL + 1 : Mesh.numX, :)) / (Mesh.numT - 1);
        VS = trapz(sol(V_SIndices,:)) / (Mesh.numS - 1);
        %VS = trapz(sol(V_SIndices,end)) / (Mesh.numS - 1);
        Msb = 2 * VS * (Geometry.SA/2) * Geometry.h_S * V.kB;
        %Msbend = Msb(end);
        %MsbendDiff = (Msbend / V.Vb) - (I.rho .* sol(I_BIndices, end));

        %VS = trapz(Mesh.x(V_TIndices), sol(V_TIndices,:));
        %Msb = 2 * VS * Geometry.W * Geometry.L * V.kB;

        dsoldt(ii,:) = (Msb ./ (V.Vb)) ...
                     + I.rho .* sol(I_BIndices, :) ... % production by I
                     - V.kL .* sol(ii, :); % loss in blood
                     
    
    %%% T in blood
    ii = T_BIndices;
                    
        dsoldt(ii,:) = T.lambda_B - T.dT_B .* sol(ii,:) ...  
                       - I.w .* sol(I_BIndices) .* sol(ii,:)...
                       - V.beta .* sol(V_BIndices,:) .* sol(ii,:);
                       ...
                       %Uavg .* sol(ii,:) .* (( V.Vb / Vb_tissue ) .* T.Ku);
                       %V.beta .* sol(V_BIndices,:) .* sol47
                       % (ii,:);
                       %Uavg .* sol(ii,:) .* T.Ku;
                       %V.beta .* sol(V_BIndices,:) .* sol(ii,:);
                       
    
    %%% I in blood
    
    ii = I_BIndices;
    
        %Iavg = trapz(sol(I_TIndices, :)) / (Mesh.numT - 1);
        %Msb = 2 * Iavg * Geometry.W * Geometry.L * Geometry.h_T * I.kI_TB;
    
        dsoldt(ii,:) =  V.beta .* sol(V_BIndices,:) .* sol(T_BIndices,:) ...
                        + I.w .* sol(ii,:) .* sol(T_BIndices,:)... 
                        - I.dI_B .* sol(ii,:); ...
                        %((1/V.Vb) .* Msb);
                        %I.dI_B .* sol(ii,:);
    
    if (includeOutputs == 1)                    
    
        cellcellInfections = I.w .* sol(I_SIndices,:) .* sol(T_SIndices,:);
        viruscellInfections = V.beta .* sol(V_SIndices, :) .* sol(T_SIndices,:);
    
        ii = VirusCellStroma;
            dsoldt(ii,:) = trapz(viruscellInfections)/(Mesh.numS - 1);
        
        ii = CellCellStroma;
            dsoldt(ii,:) = trapz(cellcellInfections)/(Mesh.numS - 1);
        
        ii = VirusCellBlood;
            dsoldt(ii,:) = I.w .* sol(I_BIndices) .* sol(T_BIndices,:);
        
        ii = CellCellBlood;
            dsoldt(ii,:) = V.beta .* sol(V_BIndices,:) .* sol(T_BIndices,:);

    end
        
    DVDT = dsoldt;
    
end 