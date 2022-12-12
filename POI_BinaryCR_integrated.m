
function [isInfected] = POI_BinaryCR_integrated(params, makePlots)

%% PARAMETERS, GEOMETRY & INITIAL CONDITIONS
% V - virus, 
% T - target cell, 
% I - infected cell, 
% U - chemokine
[V, T, I, U, Geometry] = initializeParameters(params);

% Infectivity beta model in tissue:
%getPhysics = @getPhysicsLib.physicsDiffusion_tissueInfectivity_v1;
% Collision theory model in tissue:
%getPhysics = @getPhysicsLib.physicsDiffusion_tissueCollision_v1;
% COLLISION MODEL TISSUE WITHOUT EPITHELIUM
%getPhysics = @getPhysicsLib.physicsDiffusion_tissueCollision_v1;
%getMesh = @getMeshLib.createMesh_LT;
%getIC = @getICLib.getIC_tissue;
%getResults = @getDestructuredResults.getResults_LT;
% COLLISION MODEL TISSUE WITH EPITHELIUM
getPhysics = @getPhysicsLib.physicsDiffusion_tissueCollision_v3;
getMesh = @getMeshLib.createMesh_LES;
getIC = @getICLib.getIC_EpiStroma;
getResults = @getDestructuredResults.getResults_LES;
getSparsity = @getSparsityLib.defaultSparsity;
%getSparsity = @getSparsityLib.noDrugCellMigration;
getJacobian = @getJacobianLib.jacobianNoDrugCellMigration;


%% GEOMETRY SPATIAL DISCRETIZATION

meshSize = 800;
[Mesh, totalSize] = getMesh(Geometry, meshSize);
%totalSize = Mesh.numX + Mesh.numT + Mesh.numT + 3; 
% WITHOUT CHEMOKINE:
%totalSize = Mesh.numX + Mesh.numS + Mesh.numS + 3; 
%tissue: virus: numx, chemokine: numt tcell: numt, icell: numt 
%blood: virus: 1, tcell: 1, icell: 1
     
%% TEMPORAL PARAMETRIZATION

numTimeSteps = 1000;
finalT = 20*24*60*60; %3 weeks in seconds
tvec = linspace(0,finalT,numTimeSteps);
tspan = [min(tvec),max(tvec)];
% t_sample = [0:1:10*60,...
%             10*60+1:60:24*60*60,...
%             24*60*60+1:1*60*60:36*60*60,...
%             36*60*60+1:4*60*60:tspan(end)];

%% PHYSICS - APPLYING INITIAL CONDITIONS TO MESH

IC = getIC(Mesh, V, T);
  
if (totalSize ~= length(IC))
    fprintf('Size mismatch between mesh and initial condition. \n');
    fprintf('size mesh: %d \nSize IC: %d\n', totalSize, length(IC));
    return;
end
  
%% OPTIMIZING SOLVER

S = getSparsity(Mesh, totalSize);
S = sparse(S);
%opts1 = odeset('JPattern',S);
%'RelTol',1e-3,'AbsTol',1e-4
opts1 = odeset('Vectorized','on','JPattern',S);
%opts1 = odeset('Vectorized','on','Jacobian', @(t,sol) getJacobian(t,sol,Geometry,Mesh,V,T,I,totalSize));

%% RUNNING SOLVER
      
[t,sol] = ode15s(@(t,sol) getPhysics(t,sol,IC,Geometry,Mesh,V,T,I,U), tspan, IC, opts1); %Geometry
sol = sol';

%Virions & Cells
%tissue: virus: numx, chemokine: numt tcell: numt, icell: numt 
%blood: virus: 1, tcell: 1, icell: 1

%% DECOMPACTING SOLUTION

% vir = sol(1 : Mesh.numX, :);
% chem = sol(Mesh.numX + 1 : Mesh.numX + Mesh.numT, :);
% tar = sol(Mesh.numX + Mesh.numT + 1 : Mesh.numX + (2 * Mesh.numT), :);
% inf = sol(Mesh.numX + (2 * Mesh.numT) + 1 : Mesh.numX + (3 * Mesh.numT), :);
% 
% vir_b = sol(Mesh.numX + 3 * Mesh.numT + 1, :);
% tar_b = sol(Mesh.numX + 3 * Mesh.numT + 2, :);
% inf_b = sol(Mesh.numX + 3 * Mesh.numT + 3, :);
% 
% Uavg = trapz(chem) / (Mesh.numT - 1);
% Vb_tissue = 2 * Geometry.W * Geometry.L * Geometry.h_T;
% sigmapr_term = Uavg .* tar_b .* ((Vb_tissue / V.Vb) .* T.Ku);
% r = (chem .^ 2) ./ (1 + chem + (chem .^ 2));
% r_term = U.kT .* r .* tar;
% r_avg = trapz(r_term) / (Mesh.numT - 1);

%Uavg = trapz(chem) / (Mesh.numT - 1);
%Vb_tissue = 2 * Geometry.W * Geometry.L * Geometry.h_T;
% sigmapr_term = Uavg .* tar_b .* ((Vb_tissue / V.Vb) .* T.Ku);
% r = (chem .^ 2) ./ (1 + chem + (chem .^ 2));
% r_term = U.kT .* r .* tar;
% r_avg = trapz(r_term) / (Mesh.numT - 1);

%dudx = gradient(chem) / Geometry.h_T;
%d2udx2 = gradient(dudx) / Geometry.h_T;
%dtdx = gradient(tar) / Geometry.h_T;

%stax_term = T.Stax .* T.Sigma .* (d2udx2 .* tar + dudx .* dtdx); %chemotaxis 2nd deriv
%stax_term = T.Stax .* T.Sigma .* (dudx .* tar + chem .* dtdx); %chemotaxis 1st deriv
%stax_avg = trapz(stax_term) / (Mesh.numT - 1);

[vir, tar, inf, chem, vir_b, tar_b, inf_b] = getResults(sol, Mesh);

if makePlots == 1
    settingLabels = ["showVvsTime", "showVvsX", "showUvsTime", "showCellsvsTime","showVvsTime_b","showCellsvsTime_b"];
    settingValues = [1, 1, 0, 1, 1, 1];
    settings = dictionary(settingLabels, settingValues);
    createLinePlots(t/3600, Mesh, T, vir, chem, tar, inf, vir_b, tar_b, inf_b, settings)
    %checkTermMagnitudes(t/3600, sigmapr_term, r_avg, stax_avg)

end

timeInd2Weeks = find(t >= 14*24*60*60,1);
if ~isempty(timeInd2Weeks)
    isInfected = vir_b(timeInd2Weeks) > 20; %limits of detection and limits of infection, weld et al, 2022
else
    isInfected = -1;
end



end

% function q = calcQ(cTFVDP)
% 
% IC50 = 180; %(th˚is is in ng/ml, roughly = 0.4 uM, given molar mass of 447.173 g/mol))
% %     if (IC(numx+1) ~= 0)
% q = 1./(1+(IC50./(cTFVDP)));
% 
% end

% function S = calcsparsity(mesh, totalSize)
% B = ones(totalSize,totalSize);
% % B(numx+1:2*numx+indE+4*indS,1:indG+indE) = 0; %a
% % B(1:indG+indE,numx+1:2*numx+indE+4*indS) = 0; %a
% % B(numx+1:numx+indG,2*numx+1:2*numx+indE+4*indS) = 0; %b
% % B(2*numx+1:2*numx+indE+4*indS,numx+1:numx+indG) = 0; %b
% % B(2*numx+indE+1:2*numx+indE+4*indS,numx+indG+1:numx+indG+indE) = 0; %c
% % B(numx+indG+1:numx+indG+indE,2*numx+indE+1:2*numx+indE+4*indS) = 0; %c
% % B(2*numx+indE+indS+1:2*numx+indE+4*indS,numx+indG+1:numx+indG+indE)=0; %d
% % B(numx+indG+1:numx+indG+indE,2*numx+indE+indS+1:2*numx+indE+4*indS)=0; %d
% % B(2*numx+indE+1:2*numx+indE+4*indS,2*numx+1:2*numx+indE)=0; %e
% % B(2*numx+1:2*numx+indE,2*numx+indE+1:2*numx+indE+4*indS)=0; %e
% % B(2*numx+1:2*numx+indE,numx+indG+indE+1:2*numx)=0; %f
% % B(numx+indG+indE+1:2*numx,2*numx+1:2*numx+indE)=0; %f
% 
% S = B;
% 
% end
