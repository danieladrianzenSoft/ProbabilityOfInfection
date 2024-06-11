
function [isInfected,varargout] = POI_BinaryCR_integrated(params, makePlots, includeOutputs)
    arguments
       params,
       makePlots,
       includeOutputs {mustBeNumeric} = 0
    end

    %% PARAMETERS, GEOMETRY & INITIAL CONDITIONS
    % V - virus, 
    % T - target cell, 
    % I - infected cell, 
    % U - chemokine
    % D - drug
    [V, T, I, U, D, Geometry] = initializeParameters(params);
    drugApplied = D.C0 ~= 0;
    %modifiedDrug = 0;
    
    % Infectivity beta model in tissue:
    %getPhysics = @getPhysicsLib.physicsDiffusion_tissueInfectivity_v1;
    % Collision theory model in tissue:
    %getPhysics = @getPhysicsLib.physicsDiffusion_tissueCollision_v1;
    % COLLISION MODEL TISSUE WITHOUT EPITHELIUM
    %getPhysics = @getPhysicsLib.physicsDiffusion_tissueCollision_v1;
    %getMesh = @getMeshLib.createMesh_LT;
    %getIC = @getICLib.getIC_tissue;
    %getResults = @getDestructuredResults.getResults_LT;
    
    if drugApplied == 0
        % COLLISION MODEL TISSUE WITH EPITHELIUM, NO DRUG     
        getPhysics = @getPhysicsLib.physicsDiffusion_NoDrug;
        getMesh = @getMeshLib.createMesh_LES;
        getIC = @getICLib.getIC_EpiStroma;
        getResults = @getDestructuredResults.getResults_LES;
        getSparsity = @getSparsityLib.defaultSparsity;
        %getSparsity = @getSparsityLib.noDrugCellMigration;
        %getJacobian = @getJacobianLib.jacobianNoDrugCellMigration;
    else
        % COLLISION MODEL TISSUE WITH EPITHELIUM, WITH DRUG
        %if modifiedDrug == 1
            % Drug reduces rate of cell-to-cell infection
        %    getPhysics = @getPhysicsLib.physicsDiffusion_tissueCollision_v3_Drug_cellToCell;
        %else
            % Drug reduces rate of virus-to-cell infection
            %  getPhysics = @getPhysicsLib.physicsDiffusion_tissueCollision_v3_Drug_freeVirus;
        %end
            %getPhysics = @getPhysicsLib.physicsDiffusion_tissueCollision_v3_Drug_freeVirus_wOutputs;
        switch (params.drugType)
            case 1
                getPhysics = @getPhysicsLib.physicsDiffusion_Gel;
                getMesh = @getMeshLib.createMesh_LES_Gel;
                getIC = @getICLib.getIC_EpiStroma_Gel;
                getResults = @getDestructuredResults.getResults_LES_Gel;
                getSparsity = @getSparsityLib.defaultSparsity;
            case 2
                getPhysics = @getPhysicsLib.physicsDiffusion_IVR;
                getMesh = @getMeshLib.createMesh_LES_IVR;
                getIC = @getICLib.getIC_EpiStroma_IVR;
                getResults = @getDestructuredResults.getResults_LES_IVR;
                getSparsity = @getSparsityLib.defaultSparsity;
        end
    end
    
    
    %% GEOMETRY SPATIAL DISCRETIZATION
    
    meshSize = 800;
    %[Mesh, totalSize] = getMesh(Geometry, meshSize);
    %totalSize = Mesh.numX + Mesh.numT + Mesh.numT + 3; 
    % WITHOUT CHEMOKINE:
    %totalSize = Mesh.numX + Mesh.numS + Mesh.numS + 3; 
    %tissue: virus: numx, chemokine: numt tcell: numt, icell: numt 
    %blood: virus: 1, tcell: 1, icell: 1
         
    %% TEMPORAL PARAMETRIZATION
    
    numTimeSteps = 1000;
    finalT = 20*24*60*60; %20 days in seconds
    tvec = linspace(0,finalT,numTimeSteps);
    tspan = [min(tvec),max(tvec)];
    % t_sample = [0:1:10*60,...
    %             10*60+1:60:24*60*60,...
    %             24*60*60+1:1*60*60:36*60*60,...
    %             36*60*60+1:4*60*60:tspan(end)];
    
      
    %% OPTIMIZING SOLVER
    
    % S = getSparsity(Mesh, totalSize);
    % S = sparse(S);
    % opts1 = odeset('Vectorized','on','JPattern',S,'NonNegative',true);
    %opts1 = odeset('JPattern',S);
    %'RelTol',1e-3,'AbsTol',1e-4
    %opts1 = odeset('Vectorized','on','JPattern',S,'RelTol',1e-3,'AbsTol',1e-4,'NonNegative',true);
    %opts1 = odeset('Vectorized','on','NonNegative',true);
    %opts1 = odeset('Vectorized','on','Jacobian', @(t,sol) getJacobian(t,sol,Geometry,Mesh,V,T,I,totalSize));
    
    %% RUNNING SOLVER
          
    if params.T_VD == 0
    
        % update parameters in case changes are needed due to presence or
        % absence of drug
        [V,T,I,U,D,Geometry] = updateParameters(tspan,V,T,I,U,D,Geometry,params);
        % get mesh based on updated geometry
        [Mesh, totalSize] = getMesh(Geometry, meshSize, includeOutputs);
        % get ic based on updated mesh
        IC = getIC(Mesh,V,T,D,params.T_VD, includeOutputs);
       
        if (totalSize ~= length(IC))
            fprintf('Size mismatch between mesh and initial condition. \n');
            fprintf('size mesh: %d \nSize IC: %d\n', totalSize, length(IC));
            return;
        end
        % update sparsity pattern based on updated mesh
        opts1 = odeset('Vectorized','on','JPattern',sparse(getSparsity(Mesh,totalSize)),'NonNegative',true);
    
        [t,sol] = ode15s(@(t,sol) getPhysics(t,sol,IC,Geometry,Mesh,V,T,I,D,U,includeOutputs), tspan, IC, opts1); %Geometry
        sol(sol<0) = 0;
    elseif params.T_VD < 0 || params.T_VD > 0
        % get relevant tspans depending on T_VD
        tspan_pre = [0,abs(params.T_VD)];
        tspan_post = [abs(params.T_VD),finalT];
    
        % get parameters at phase 1 or 2 of simulation (depending on T_VD)
        [V_pre,T_pre,I_pre,U_pre,D_pre,Geometry_pre] = updateParameters(tspan_pre,V,T,I,U,D,Geometry,params);
        [V_post,T_post,I_post,U_post,D_post,Geometry_post] = updateParameters(tspan_post,V,T,I,U,D,Geometry,params);
    
        % construct mesh based on geometry with highest h_L to maximize 
        % resolution in this region
        if (Geometry_pre.h_L > Geometry_post.h_L)
            [Mesh, totalSize] = getMesh(Geometry_pre, meshSize, includeOutputs);
        else
            [Mesh, totalSize] = getMesh(Geometry_post, meshSize, includeOutputs);
        end
    
        % update sparsity pattern based on constructed mesh
        opts1 = odeset('Vectorized','on','JPattern',sparse(getSparsity(Mesh,totalSize)),'NonNegative',true);
    
        % get ic, given chosen mesh, at both phases of simulation
        IC_pre = getIC(Mesh,V,T,D,params.T_VD, includeOutputs);
    
        if (totalSize ~= length(IC_pre))
            fprintf('Size mismatch between mesh and initial condition. \n');
            fprintf('Size mesh: %d \nSize IC: %d\n', totalSize, length(IC_pre));
            return;
        end
    
        % run solution for 1st phase
        [t1,sol1] = ode15s(@(t,sol) getPhysics(t,sol,IC_pre,Geometry_post,Mesh,V_pre,T_pre,I_pre,D_pre,U_pre, includeOutputs), tspan_pre, IC_pre, opts1); %Geometry
        
        % get ic for 2nd phase of simulation based on prior solution
        switch (params.drugType)
            case 1
                IC_post = getICLib.getIC_EpiStroma_Gel_TVD(Mesh,V_post,T_post,D_post,params.T_VD,sol1,getResults, includeOutputs);
            case 2
                IC_post = getICLib.getIC_EpiStroma_IVR_TVD(Mesh,V_post,T_post,D_post,params.T_VD,sol1,getResults, includeOutputs);
        end
    
        % run solution for 2nd phase
        [t2,sol2] = ode15s(@(t,sol) getPhysics(t,sol,IC_post,Geometry_post,Mesh,V_post,T_post,I_post,D_post,U_post, includeOutputs), tspan_post, IC_post, opts1); %Geometry
    
        t = [t1(1:end-1);t2];
        sol = [sol1(1:end-1,:);sol2];  
        sol(sol<0) = 0;
    else
        fprintf('\nIssue with T_VD. \n');
    end
    sol = sol';
    
    %Virions & Cells
    %tissue: virus: numx, chemokine: numt tcell: numt, icell: numt 
    %blood: virus: 1, tcell: 1, icell: 1
    
    %% DECOMPACTING SOLUTION
    
    if includeOutputs == 0
        [vir, drug, drugdp, tar, inf, chem, vir_b, drug_b, drugdp_b, tar_b, inf_b, ~] = getResults(sol, Mesh, 1, includeOutputs);
    else
        [vir, drug, drugdp, tar, inf, chem, vir_b, drug_b, drugdp_b, tar_b, inf_b, vir_tar_t, inf_tar_t, vir_tar_b, inf_tar_b] = getResults(sol, Mesh, 1, includeOutputs);
         tInterp = linspace(0,finalT,40);
         vir_tar_t_interp = interp1(t,vir_tar_t,tInterp,'spline');
         inf_tar_t_interp = interp1(t,inf_tar_t,tInterp,'spline');
         vir_tar_b_interp = interp1(t,vir_tar_b,tInterp,'spline');
         inf_tar_b_interp = interp1(t,inf_tar_b,tInterp,'spline');
         
         varargout{1} = tInterp;
         varargout{2} = vir_tar_t_interp;
         varargout{3} = inf_tar_t_interp;
         varargout{4} = vir_tar_b_interp;
         varargout{5} = inf_tar_b_interp;
    end
    
    if makePlots == 1
        settingLabels = ["showVvsTime", "showVvsX", "showDvsTime", "showDdpvsTime", "showUvsTime", "showCellsvsTime","showVvsTime_b", "showDvsTime_b","showDdpvsTime_b","showqvsTime_b","showCellsvsTime_b"];
        if drugApplied == 0 
            settingValues = [1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 1];
        else
            settingValues = [1, 1, 1, 1, 0, 1, 1, 1, 1, 0, 1];
        end
        %settingValues = [0, 0, 1, 1, 0, 0, 0, 1, 1, 1, 0];
        settings = dictionary(settingLabels, settingValues);
        createLinePlots(t/3600, Mesh, T, vir, drug, drugdp, chem, tar, inf, vir_b, drug_b, drugdp_b, tar_b, inf_b, params.drugType, settings)

        if (includeOutputs ~= 0)
            additionalOutputSettingLabels = ["showCellCellAndVirusCellAllCompartments","showCellCellOverVirusCellPerCompartment", "showInfsInTissueOverBlood"];
            additionalOutputSettingValues = [1,1,1];
            additionalOutputSettings = dictionary(additionalOutputSettingLabels, additionalOutputSettingValues);

            createAdditionalOutputPlots(t/3600, vir_tar_t, inf_tar_t, vir_tar_b, inf_tar_b, additionalOutputSettings)
        end
        %createHeatmap(t/3600, Mesh, drug, [zeros(Mesh.numL,length(t));drugdp], vir, drugApplied)
        %checkTermMagnitudes(t/3600, sigmapr_term, r_avg, stax_avg)
    
    end
    
    timeInd2Weeks = find(t >= 18*24*60*60,1);
    if ~isempty(timeInd2Weeks)
        isInfected = vir_b(timeInd2Weeks) > 20; %limits of detection and limits of infection, weld et al, 2022
        % 1 if > 20, 0 if < 20
    else
        isInfected = -1;
        % error, no time point greater than 2 weeks.
    end



end

% function q = calcQ(cTFVDP)
% 
% IC50 = 180; %(th˚is is in ng/ml, roughly = 0.4 uM, given molar mass of 447.173 g/mol))
% %     if (IC(numx+1) ~= 0)
% q = 1./(1+(IC50./(cTFVDP)));
% 
% end