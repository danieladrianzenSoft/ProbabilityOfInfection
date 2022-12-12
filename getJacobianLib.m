classdef getJacobianLib
    methods(Static)
        function J = jacobianNoDrugCellMigration(t,sol,Geometry,Mesh,V,T,I,totalSize)
            J = zeros(totalSize,totalSize);
            dxL = Geometry.h_L / (Mesh.numL - 1);
            dxE = Geometry.h_E / (Mesh.numE - 1);
            dxS = Geometry.h_S / (Mesh.numS - 1);

            V_SIndices = Mesh.numL + Mesh.numE + 1 : Mesh.numX;
            T_SIndices = Mesh.numX + 1 : Mesh.numX + Mesh.numS;
            I_SIndices = Mesh.numX + Mesh.numS + 1 : Mesh.numX + 2 * Mesh.numS;
            V_BIndices = Mesh.numX + 2 * Mesh.numS + 1;
            T_BIndices = Mesh.numX + 2 * Mesh.numS + 2;
            I_BIndices = Mesh.numX + 2 * Mesh.numS + 3;

            for i = 1:totalSize
                if i == 1
                    J(i,i) = -2*V.DV_L / (dxL^2) - V.kD;
                    J(i,i+1) = 2*V.DV_L / (dxL^2);
                elseif (i > 1 && i < Mesh.numL)
                    J(i,i-1) = V.DV_L / (dxL^2);
                    J(i,i) = -2*V.DV_L / (dxL^2) - V.kD;
                    J(i,i+1) = V.DV_L / (dxL^2);
                elseif i == Mesh.numL
                    J(i,i-1) = V.DV_L / (dxL^2);
                    J(i,i) = (V.DV_L / (dxL^2)) * (-2 + V.DV_L/(V.DV_E*V.phi_LT+V.DV_L)) - V.kD;
                    J(i,i+1) = (V.DV_L / (dxL^2)) * (V.DV_E/(V.DV_E*V.phi_LT+V.DV_L));
                elseif i == Mesh.numL+1
                    J(i,i-1) = V.DV_E / (dxE^2) * (V.DV_L / (V.DV_L/V.phi_LT + V.DV_E));
                    J(i,i) = V.DV_E / (dxE^2) * ((V.DV_E / (V.DV_L/V.phi_LT + V.DV_E))-2);
                    J(i,i+1) = V.DV_E / (dxE^2);
                elseif (i > Mesh.numL+1 && i < Mesh.numL + Mesh.numE)
                    J(i,i-1) = V.DV_E / (dxE^2);
                    J(i,i) = -2*(V.DV_E / (dxE^2));
                    J(i,i+1) = V.DV_E / (dxE^2);
                elseif i == Mesh.numL+Mesh.numE
                    J(i,i-1) = V.DV_E / (dxE^2);
                    J(i,i) = (V.DV_E / (dxE^2))*(V.DV_E/(V.DV_S+V.DV_E) - 2);
                    J(i,i+1) = (V.DV_E / (dxE^2))*(V.DV_S/(V.DV_S+V.DV_E));
                elseif i == V_SIndices(1)
                    J(i,i-1) = (V.DV_S / (dxS^2))*(V.DV_E / (V.DV_E+V.DV_S));
                    J(i,i) = (V.DV_S / (dxS^2))*(-2 + (V.DV_S / (V.DV_E+V.DV_S))) - (V.kB);
                    J(i,i+1) = (V.DV_S / (dxS^2));
                    J(i,i+2*Mesh.numS) = I.rho;
                elseif (i > V_SIndices(1) && i < V_SIndices(end))
                    J(i,i-1) = V.DV_S / (dxS^2);
                    J(i,i) = -2*(V.DV_S / (dxS^2)) - (V.kB);
                    J(i,i+1) = (V.DV_S / (dxS^2));
                    J(i,i+2*Mesh.numS) = I.rho;
                elseif i == V_SIndices(end)
                    J(i,i-1) = 2*V.DV_S / (dxS^2);
                    J(i,i) = -2*(V.DV_S / (dxS^2)) - (V.kB);
                    J(i,i+1) = 0;
                    J(i,i+2*Mesh.numS) = I.rho;
                elseif i == T_SIndices(1)
                    J(i,i-Mesh.numS) = -V.beta*sol(i,:);
                    J(i,i-1) = 0;
                    J(i,i) = -2*(T.DT_T / (dxS^2)) - T.dT_T - I.w*sol(i+Mesh.numS,:) - V.beta*sol(i-Mesh.numS,:);
                    J(i,i+1) = 2*(T.DT_T / (dxS^2));
                    J(i,i+Mesh.numS) = -I.w*sol(i,:);
                elseif (i > T_SIndices(1) && i < T_SIndices(end))
                    J(i,i-Mesh.numS) = -V.beta*sol(i,:);
                    J(i,i-1) = T.DT_T / dxS^2;
                    J(i,i) = -2*(T.DT_T / (dxS^2)) - T.dT_T - I.w*sol(i+Mesh.numS,:) - V.beta*sol(i-Mesh.numS,:);
                    J(i,i+1) = T.DT_T / (dxS^2);
                    J(i,i+Mesh.numS) = -I.w*sol(i,:);
                elseif i == T_SIndices(end)
                    J(i,i-Mesh.numS) = -V.beta*sol(i,:);
                    J(i,i-1) = 2*T.DT_T / dxS^2;
                    J(i,i) = -2*(T.DT_T / (dxS^2)) - T.dT_T - I.w*sol(i+Mesh.numS,:) - V.beta*sol(i-Mesh.numS,:);
                    J(i,i+1) = 0;
                    J(i,i+Mesh.numS) = -I.w*sol(i,:);
                elseif i == I_SIndices(1)
                    J(i,i-2*Mesh.numS) = V.beta*sol(i-Mesh.numS,:);
                    J(i,i-Mesh.numS) = V.beta*sol(i-2*Mesh.numS,:)+I.w*sol(i,:);
                    J(i,i-1) = 0;
                    J(i,i) = -2*(I.DI_T / (dxS^2)) + I.w*sol(i-Mesh.numS,:) - I.dI_T;
                    J(i,i+1) = 2*(I.DI_T / (dxS^2));
                elseif (i > I_SIndices(1) && i < I_SIndices(end))
                    J(i,i-2*Mesh.numS) = V.beta*sol(i-Mesh.numS,:);
                    J(i,i-Mesh.numS) = V.beta*sol(i-2*Mesh.numS,:)+I.w*sol(i,:);
                    J(i,i-1) = I.DI_T / (dxS^2);
                    J(i,i) = -2*(I.DI_T / (dxS^2)) + I.w*sol(i-Mesh.numS,:) - I.dI_T;
                    J(i,i+1) = I.DI_T / (dxS^2);
                elseif i == I_SIndices(end)
                    J(i,i-2*Mesh.numS) = V.beta*sol(i-Mesh.numS,:);
                    J(i,i-Mesh.numS) = V.beta*sol(i-2*Mesh.numS,:)+I.w*sol(i,:);
                    J(i,i-1) = 2*I.DI_T / (dxS^2);
                    J(i,i) = -2*(I.DI_T / (dxS^2)) + I.w*sol(i-Mesh.numS,:) - I.dI_T;
                    J(i,i+1) = 0;
                elseif i == V_BIndices
                    %J(i,i-3*Mesh.numS:i-2*Mesh.numS-1) = 2*(Geometry.W*Geometry.L / V.Vb) * V.kB * sol(i-3*Mesh.numS:i-2*Mesh.numS-1,:);
                    J(i,i-3*Mesh.numS:i-2*Mesh.numS-1) = (2*V.kB*Geometry.W*Geometry.L*Geometry.h_S) / (V.Vb * Mesh.numS);
                    J(i,i) = -V.kL;
                    J(i,i+1) = 0;
                    J(i,i+2) = I.rho;
                elseif i == T_BIndices
                    J(i,i-1) = -V.beta*sol(i,:);
                    J(i,i) = -T.dT_B - I.w*sol(i+1,:) - V.beta*sol(i-1,:);
                    J(i,i+1) = -I.w*sol(i);
                elseif i == I_BIndices
                    J(i,i-2) = V.beta*sol(i-1,:);
                    J(i,i-1) = V.beta*sol(i-2,:) + I.w*sol(i,:);
                    J(i,i) = I.w*sol(i-1,:) - I.dI_B;
                end
            end
        end
        function DVDT = jacSymbolicNoDrugCellMigration(IC,Geometry,Mesh,V,T,I)
            syms t
            sol = sym('sol', [1 length(IC)]);
            dsoldt = sym('dsoldt', [1 length(IC)]);
            dsoldx = sym('dsoldx', [1 length(IC)]);
            d2soldx2 = sym('d2soldx2', [1 length(IC)]);

%             sym ('dcdt', []);
%             dsoldt = zeros(length(IC),1);
%             dsoldx = zeros(length(IC),1);
%             d2soldx2 = zeros(length(IC),1);
            
            %zeta = zeros(Mesh.indT, size(sol,2));
            dxL = Geometry.h_L / (Mesh.numL - 1);
            dxE = Geometry.h_E / (Mesh.numE - 1);
            dxS = Geometry.h_S / (Mesh.numS - 1);

            %dxL = Mesh.x(Mesh.numL)-Mesh.x(Mesh.numL-1);
            %dxT = Mesh.x(Mesh.numT)-Mesh.x(Mesh.numT-1);

            %HIV INTERFACE
            %V_intf1a = (V.DV_L.*sol(Mesh.numL)+(V.DV_T.*sol(Mesh.numL+1)))./((V.DV_T.*V.phi_LT)+V.DV_L); %V at lumen/tissue interface
            %V_intf1b = (V.DV_L.*sol(Mesh.numL)+(V.DV_T.*sol(Mesh.numL+1)))./((V.DV_L./V.phi_LT)+V.DV_T);
            V_intf1a = (V.DV_L.*sol(Mesh.numL)+(V.DV_E.*sol(Mesh.numL+1)))./((V.DV_E.*V.phi_LT)+V.DV_L); %V at lumen/epithelium interface
            V_intf1b = (V.DV_L.*sol(Mesh.numL)+(V.DV_E.*sol(Mesh.numL+1)))./((V.DV_L./V.phi_LT)+V.DV_E);
            V_intf2a = (V.DV_E.*sol(Mesh.numL+Mesh.numE)+(V.DV_S.*sol(Mesh.numL+Mesh.numE+1)))./(V.DV_S+V.DV_E); % V at epithelium/tissue interface, assuming part.coeff of 1
            V_intf2b = (V.DV_E.*sol(Mesh.numL+Mesh.numE)+(V.DV_S.*sol(Mesh.numL+Mesh.numE+1)))./(V.DV_E+V.DV_S);
            
            %INDICES
            V_EIndices = Mesh.numL + 1 : Mesh.numL + Mesh.numE;
            V_SIndices = Mesh.numL + Mesh.numE + 1 : Mesh.numX;
            T_SIndices = Mesh.numX + 1 : Mesh.numX + Mesh.numS;
            I_SIndices = Mesh.numX + Mesh.numS + 1 : Mesh.numX + 2 * Mesh.numS;
            V_BIndices = Mesh.numX + 2 * Mesh.numS + 1;
            T_BIndices = Mesh.numX + 2 * Mesh.numS + 2;
            I_BIndices = Mesh.numX + 2 * Mesh.numS + 3;

            %T_TIndices = Mesh.numX + 1 : (Mesh.numX + Mesh.numT);
            %I_TIndices = (Mesh.numX + Mesh.numT + 1) : Mesh.numX + (2 * Mesh.numT); 
            % U_BIndices = Mesh.numX + (3 * Mesh.numT) + 1;
            %V_BIndices = Mesh.numX + (2 * Mesh.numT) + 1;
            %T_BIndices = Mesh.numX + (2 * Mesh.numT) + 2;
            %I_BIndices = Mesh.numX + (2 * Mesh.numT) + 3;
            

            %% %%%%%%%%% TISSUE %%%%%%%%%% %%
            
            %%%%% VIRUS %%%%%%
            
            ii = 1; %y direction BC zero flux (dcdx=0)
                d2soldx2(ii) = 2 * (sol(ii+1) - sol(ii))./(dxL.^2);
                dsoldt(ii) = V.DV_L .* d2soldx2(ii) - (V.kD.*sol(ii));
                %dVdt(i) = 2*D_G.*(V(i+1)-V(i))./(dx.^2);
                
            ii = 2:Mesh.numL-1; %in lumen
                %dVdx(ii) = (V(ii+1)-V(ii))./dxG;
                d2soldx2(ii) = (sol(ii+1) - 2.*sol(ii) + sol(ii-1))./(dxL.^2);
                dsoldt(ii) = V.DV_L.*d2soldx2(ii) - (V.kD.*sol(ii));
            
            ii = Mesh.numL; %right before interface - lumen
                %dVdx(ii) = (V_intf1a - V(ii))./dxG;
                d2soldx2(ii) = (V_intf1a - 2.*sol(ii) + sol(ii-1))./(dxL.^2);
                dsoldt(ii) = V.DV_L.*d2soldx2(ii) - (V.kD.*sol(ii));
            
            ii = V_EIndices(1); %right after interface - in epithelium
                %dVdx(ii) = (V(ii+1)-V(ii))./dxE;
                d2soldx2(ii) = (sol(ii+1) - 2.*sol(ii) + (V_intf1b))./(dxE.^2);
                dsoldt(ii) = V.DV_E.*d2soldx2(ii) - (V.kL).*sol(ii);
            
            ii = V_EIndices(2:end-1); %in epithelium
                %dVdx(ii) = (V(ii+1)-V(ii))./dxE;
                d2soldx2(ii) = (sol(ii+1) - 2.*sol(ii) + sol(ii-1))./(dxE.^2);
                dsoldt(ii) = V.DV_E.*d2soldx2(ii) - (V.kL).*sol(ii);
            
            ii = V_EIndices(end); %right before interface - in epithelium
                %dVdx(ii) = (V_intf2a-V(ii))./dxE;
                d2soldx2(ii) = (V_intf2a - 2.*sol(ii) + sol(ii-1))./(dxE.^2);
                dsoldt(ii) = V.DV_E*d2soldx2(ii) - (V.kL).*sol(ii); %BC zero flux

            ii = V_SIndices(1); %right after interface - in stroma
                %dVdx(ii) = (V(ii+1)-V(ii))./dxE;
                d2soldx2(ii) = (sol(ii+1) - 2.*sol(ii) + (V_intf2b))./(dxS.^2);
                dsoldt(ii) = V.DV_S.*d2soldx2(ii) - (V.kB + V.kL).*sol(ii) + I.rho*sol(I_SIndices(1));
            
            ii = V_SIndices(2:end-1); %in stroma
                %dVdx(ii) = (V(ii+1)-V(ii))./dxE;
                d2soldx2(ii) = (sol(ii+1) - 2.*sol(ii) + sol(ii-1))./(dxS.^2);
                dsoldt(ii) = V.DV_S.*d2soldx2(ii) - (V.kB + V.kL).*sol(ii) + I.rho*sol(I_SIndices(2:end-1));
            
            ii = V_SIndices(end); %end of stroma
                %dVdx(ii) = (V_intf2a-V(ii))./dxE;
                d2soldx2(ii) = 2*(sol(ii-1)-sol(ii))./(dxS.^2);
                dsoldt(ii) = V.DV_S*d2soldx2(ii) - (V.kB + V.kL).*sol(ii) + I.rho*sol(I_SIndices(end)); %BC zero flux            

            %%%%% CELLS %%%%%%
            
            zeta_T = 4 * pi * (T.r + V.r) .* (T.DT_T + V.DV_T); 
        
            
            ii = T_SIndices(1); %target cells
                
                d2soldx2(ii) = 2.*(sol(ii+1)-sol(ii))./(dxS.^2);
                dsoldt(ii) = ...
                             T.DT_T.*d2soldx2(ii) ... %diffusion
                             + (T.lambda_T - T.dT_T .* sol(ii)) ... % proliferation rate
                             - I.w .* sol(I_SIndices(1)) .* sol(ii)... % loss due to cell-cell infections
                             - V.s .* zeta_T .* sol(V_SIndices(1)) .* sol(ii); % loss due to infections, collision theory
                                        
            ii = T_SIndices(2:end-1); %target cells
            
            
                dsoldx(ii) = (sol(ii+1)-sol(ii-1))./(2*dxS);
                d2soldx2(ii) = (sol(ii+1)-2.*sol(ii)+sol(ii-1))./(dxS.^2);
                dsoldt(ii) = ... 
                             T.DT_T.*d2soldx2(ii) ... %diffusion
                             + (T.lambda_T - T.dT_T .* sol(ii)) ... % proliferation rate,
                             - I.w .* sol(I_SIndices(2:end-1)) .* sol(ii)... % loss due to cell-cell infections
                             - V.s .* zeta_T .* sol(V_SIndices(2:end-1)).* sol(ii); % loss due to infections, collision theory
                             ...V.beta .* sol(V_SIndices(2:end-1), :) .* sol(ii);

                            ...T.Stax .* T.Sigma .* (dsoldx(U_TIndices(2:end-1), :) .* sol(ii, :) + dsoldx(ii) .* sol(U_TIndices(2:end-1), :)); % chemotaxis first derivative u
                             
                            %Uavg .* sol(T_BIndices) .* ((Vb_tissue / V.Vb) .* T.Ku) + ... % uptake from blood due to U(x,t)
                             % T.Stax .* T.Sigma .* (d2soldx2(U_TIndices(2:end-1)) .* sol(ii) + dsoldx(U_TIndices(2:end-1))) .* dsoldx(ii); %chemotaxis second derivative u
                             % V.s .* zeta_T .* sol((Mesh.numL+2):(Mesh.numL+Mesh.numT)-1) .* sol(ii); % loss due to infections
                             %V.s .* zeta_T .* sol((Mesh.numL+2):(Mesh.numL+Mesh.numT)-1) .* sol(ii); % loss due to infections
                             %Uavg .* sol(T_BIndices) .* T.Ku; %uptake from blood due to U(x,t)
                
            ii = T_SIndices(end); % target cells
            
                d2soldx2(ii) = 2*(sol(ii-1)-sol(ii))./(dxS.^2);
                dsoldt(ii) = ...
                             T.DT_T.*d2soldx2(ii) ... %diffusion
                             + (T.lambda_T - T.dT_T .* sol(ii)) ... % proliferation rate
                             - I.w .* sol(I_SIndices(end)) .* sol(ii) ... % loss due to cell-cell infections
                             - V.s .* zeta_T .* sol(V_SIndices(end)).* sol(ii); % loss due to infections, collision theory
            
            ii = I_SIndices(1); %infected cells
            
                d2soldx2(ii) = 2.*(sol(ii+1)-sol(ii))./(dxS.^2);
                dsoldt(ii) = ...
                             I.DI_T .* d2soldx2(ii) ... %diffusion
                             + V.s.*zeta_T .* sol(V_SIndices(1)) .* sol(T_SIndices(1)) ... %increase due to infections, collision theory
                             + I.w .* sol(ii) .* sol(T_SIndices(1)) ...
                             - I.dI_T * sol(ii); %death rate
                        
            ii = I_SIndices(2:end-1); %infected cells
            
                d2soldx2(ii) = (sol(ii+1)-2.*sol(ii)+sol(ii-1))./(dxS.^2);
                dsoldt(ii) = ... 
                             I.DI_T .* d2soldx2(ii) ... %diffusion
                             + V.s.*zeta_T .* sol(V_SIndices(2:end-1)) .* sol(T_SIndices(2:end-1)) - ... %increase due to infections, collision theory
                             ...V.beta .* sol(V_SIndices(2:end-1), :) .* sol(T_SIndices(2:end-1)) - ...
                             + I.w .* sol(ii) .* sol(T_SIndices(2:end-1)) ...
                             - I.dI_T * sol(ii); %death rate 
                             %I.dI_T * sol(ii) - ... %death rate
                             %I.kI_TB * sol(ii); %loss to blood.
            
                
            ii = I_SIndices(end); %infected cells
            
                d2soldx2(ii) = 2*(sol(ii-1)-sol(ii))./(dxS.^2);
                dsoldt(ii) = ...
                             I.DI_T.*d2soldx2(ii) ... %diffusion
                             + V.s.*zeta_T .* sol(V_SIndices(end)) .* sol(T_SIndices(end)) - ... %increase due to infections, collision theory
                             + I.w .* sol(ii) .* sol(T_SIndices(end)) ...
                             ...V.beta .* sol(V_SIndices(end), :) .* sol(T_SIndices(end)) - ...
                             - I.dI_T * sol(ii); %death rate
            
            
            %% %%%%%%%%% BLOOD %%%%%%%%%% %%
                   
            %%% V in blood
            ii = V_BIndices;
            
                VS = trapz(sol(V_SIndices)) / (Mesh.numS - 1);
                Msb = 2 * VS * Geometry.W * Geometry.L * Geometry.h_S * V.kB;

                dsoldt(ii) = (Msb ./ (V.Vb)) ...
                             + I.rho .* sol(I_BIndices) ... % production by I
                             - V.kL .* sol(ii); % loss in blood
                             
            
            %%% T in blood
            ii = T_BIndices;
                            
                dsoldt(ii) = T.lambda_B - T.dT_B .* sol(ii) ...  
                               - I.w .* sol(I_BIndices) .* sol(ii)...
                               - V.beta .* sol(V_BIndices) .* sol(ii);
                               
            
            %%% I in blood
            
            ii = I_BIndices;
            
                dsoldt(ii) =  V.beta .* sol(V_BIndices) .* sol(T_BIndices) ...
                                + I.w .* sol(ii) .* sol(T_BIndices)... 
                                - I.dI_B .* sol(ii); ...
            
            J = jacobian(dsoldt,sol);    
            DVDT = matlabFunction(J,'Vars',{t,[sol.']});

        end
    end
end