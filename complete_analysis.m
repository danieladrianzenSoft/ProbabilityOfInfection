    clear
    clc
    
    % simConditions:
    % 0 - overwrite entire dataset.
    % 1 - restart simulation from last saved batch. 
    % 2 - add new parameters to existing simulation, evaluate batches and append.
    % 3 - plot params and evaluate representative result for prior sim where isInfected = 1.
    % 4 - run independent sim with specified parameters.
    % 5 - run independent sim with parameters specified from file

    % PARAMETER DEFINITIONS:
    % phase: 0 - Midcycle, 1 - Follicular, 2 - Luteal
    % drugType: 0 - None, 1 - Gel, 2 - IVR 

    %simConditions = 0;
    %runCluster = 1;
    %readParamMeansFromFile = 1;
    %includeAdditionalOutputs = 1;

    simConditions = 0;
    runCluster = 1;
    readParamMeansFromFile = 1;
    includeAdditionalOutputs = 1;

    N = 100000;
    tictocStep = 2500;
    batchSize = 10000;

    %N = 10;
    %tictocStep = 5;
    %batchSize=5;
    
    % C0 = 0;
    % DVdelay = 0;
    %N = 200;
    %tictocStep = 50;
    %batchSize=50;
    %C0 = 0;
    %DVdelay = 0;

    %paramIndices = 1:length(paramNames);
    %paramValues = {K, r0, Rtot, Q, lambda, mu, phi, nmax};
    %params = dictionary(paramNames, paramValues);
    %paramLabels = dictionary(paramIndices,paramNames);

    dateStamp = string(datetime("now"));

    filenameParams = 'parameters_1.xlsx';
    filenameResults = 'resultsPOI_1.xlsx';
    filenameConfig = 'SimulationConfig_1.txt';
    filenameSimStatus = 'SimulationStatus.xlsx';
    filenameAdditionalOutputs = 'AdditionalOutputs_1.xlsx';

    if simConditions <= 3
        if readParamMeansFromFile == 0
            paramLabels = ["V_0", "SA","h_E","rho", "k_B", "k_L", "beta", "w", "T0_B", "T0_T", "Vf", "Vs", "Vg", "C0", "T_VD", "phase", "drugType"];
            V_0 = 1e4;
            h_E = 0.02;
            rho = 850;
            k_B = 5;
            k_L = 10;
            T0_B = 1e3;
            T0_T = 1e3;
            beta = 6.5*10^(-7);
            w = (2e-4);
            SA = 100; % surface area vaginal epithelium in humans
            Vg = 2; % mL volume of drug gel
            Vs = 2; % mL volume of semen
            Vf = 1; % mL volume of vaginal fluid
            C0 = 0;
            T_VD = 0;
            phase = 0; %PHASE = 0: MIDCYCLE, 1: FOLLICULAR, 2: LUTEAL
            drugType = 0;
            paramValues = [V_0,SA,h_E,rho,k_B,k_L,beta,w,T0_B,T0_T,Vf,Vs,Vg,C0,T_VD,phase,drugType];
            paramMeans = dictionary(paramLabels,paramValues);
            paramMeanCombinations = table(paramLabels,V_0,SA,h_E,rho,k_B,k_L,beta,w,T0_B,T0_T,Vf,Vs,Vg,C0,T_VD,phase,drugType);
        else
            paramMeanCombinations = readtable('jmp_doe.xlsx');
            paramMeanCombinations(:,1) = [];
            paramLabels = string(paramMeanCombinations.Properties.VariableNames);
        end
    end

    
    if simConditions < 3

        if simConditions == 1 || simConditions == 2
            startCombinationFile = readtable('SimulationStatus.xlsx');
            startCombination = startCombinationFile.combination;
        else
            startCombination = 1;
        end


        for combination = startCombination:size(paramMeanCombinations,1)
            if N > batchSize
                batchNum = floor(N/batchSize);
            else
                batchNum = 1;
            end

            if readParamMeansFromFile == 1
                filenameParams = strsplit(filenameParams,'.');
                filenameResults = strsplit(filenameResults,'.');
                filenameConfig = strsplit(filenameConfig,'.');
                filenameAdditionalOutputs = strsplit(filenameAdditionalOutputs,'.');
                filenameParams = strsplit(filenameParams{1},'_');
                filenameResults = strsplit(filenameResults{1},'_');
                filenameConfig = strsplit(filenameConfig{1},'_');
                filenameAdditionalOutputs = strsplit(filenameAdditionalOutputs{1},'_');

                filenameParams = sprintf("%s_%d.xlsx",filenameParams{1}, combination);
                filenameResults = sprintf("%s_%d.xlsx",filenameResults{1}, combination);
                filenameConfig = sprintf("%s_%d.txt",filenameConfig{1}, combination);
                filenameAdditionalOutputs = sprintf("%s_%d.xlsx",filenameAdditionalOutputs{1}, combination);

                paramMeans = dictionary(paramLabels,paramMeanCombinations{combination,:});
            end
                
            if simConditions == 0
                overwrite = 1;
                newparams = samplingNormalPOI_BinaryCR(filenameParams,N,paramMeans,overwrite);
                OutputPlaceholderTable = array2table(zeros(0,size(newparams,2)+1));
                OutputPlaceholderTable.Properties.VariableNames = {'V_0','SA','h_E',...
                   'rho','k_B','k_L','beta','w','T0_B','T0_T','Vf','Vs','Vg','C0','T_VD','phase','drugType','isInfected'};

                if (includeAdditionalOutputs == 1)
                   AdditionalOutputsPlaceholderTable = array2table(zeros(0,6));
                   AdditionalOutputsPlaceholderTable.Properties.VariableNames = {'sim_id','t','vir_tar_t','inf_tar_t','vir_tar_b','inf_tar_b'};
                   writetable(AdditionalOutputsPlaceholderTable,filenameAdditionalOutputs,'WriteMode','overwrite','WriteRowNames',true)
                end
                writetable(OutputPlaceholderTable,filenameResults,'WriteMode','overwritesheet','WriteRowNames',true)
                fileID = fopen(filenameConfig,'w+');
                fprintf(fileID,'\nDate Started: %s\n',dateStamp);
                fprintf(fileID,'Number of Iterations: %d\n',N);
                fprintf(fileID,'Number of Batches: %d\n',batchNum);
                startBatch = 1;
            elseif simConditions == 1
                overwrite = 0;
                try
                    newparams = readtable(filenameParams);
                    fileID = fopen(filenameConfig,'a');
                    startBatch = (size(readtable(filenameResults),1)/batchSize)+1;
                    fprintf(fileID,'\n ------------------------------- \n');
                    fprintf(fileID,'\nDate Resumed: %s\n',dateStamp);
                    fprintf(fileID,'Starting Batch: %d\n',startBatch);
                catch
                    newparams = samplingNormalPOI_BinaryCR(filenameParams,N,paramMeans,1);
                    OutputPlaceholderTable = array2table(zeros(0,size(newparams,2)+1));
                    OutputPlaceholderTable.Properties.VariableNames = {'V_0','SA','h_E',...
                       'rho','k_B','k_L','beta','w','T0_B','T0_T','Vf','Vs','Vg','C0','T_VD','phase','drugType','isInfected'};
                    if (includeAdditionalOutputs == 1)
                       AdditionalOutputsPlaceholderTable = array2table(zeros(0,6));
                       AdditionalOutputsPlaceholderTable.Properties.VariableNames = {'sim_id','t','vir_tar_t','inf_tar_t','vir_tar_b','inf_tar_b'};
                       writetable(AdditionalOutputsPlaceholderTable,filenameAdditionalOutputs,'WriteMode','overwrite','WriteRowNames',true)
                    end
                    writetable(OutputPlaceholderTable,filenameResults,'WriteMode','overwritesheet','WriteRowNames',true)
                    fileID = fopen(filenameConfig,'w+');
                    fprintf(fileID,'\nDate Started: %s\n',dateStamp);
                    fprintf(fileID,'(Resumed from prior set of simulations)\n\n');
                    fprintf(fileID,'Number of Iterations: %d\n',N);
                    fprintf(fileID,'Number of Batches: %d\n',batchNum);
                    startBatch = 1;
                end
            elseif simConditions == 2
                overwrite = 0;
                newparams = samplingNormalPOI_BinaryCR(filenameParams,N,C0,DVdelay,paramMeans,overwrite);
                completeParams = readtable(filenameParams);
                fileID = fopen(filenameConfig,'a');
                fprintf(fileID,'\n ------------------------------- \n');
                fprintf(fileID,'\nDate Restarted with New Batches: %s\n',dateStamp);
                fprintf(fileID,'Number of Additional Iterations: %d\n',N);
                fprintf(fileID,'Number of Additional Batches: %d\n',batchNum);
                fprintf(fileID,'Total Number of Iterations: %d\n',size(completeParams,1));
                startBatch = 1;        
            end
    
            cumTime = 0;
                        
            for j = startBatch:batchNum
    
                if mod(N,batchNum) ~= 0 && j == batchNum
                    numItersBatch = batchSize + mod(N,batchNum);
                else
                    numItersBatch = floor(N/batchNum);
                end
    
                start = tic;
                batchParams = newparams((j-1)*batchSize+1:j*batchSize,:);
                isInfected = zeros(numItersBatch,1);
                t_inf = zeros(numItersBatch,1);

                vir_tar_t = cell(numItersBatch,1);
                inf_tar_t = cell(numItersBatch,1);
                vir_tar_b = cell(numItersBatch,1);
                inf_tar_b = cell(numItersBatch,1);
                t_outputs = cell(numItersBatch,1);    
                sim_id = cell(numItersBatch,1);
        
                parfor i = 1:numItersBatch
                    %[isInfected(i)] = POI_BinaryCR_integrated(batchParams(i,:), 0);
                    t_outputs_batch = [];
                    vir_tar_t_batch = [];
                    inf_tar_t_batch = [];
                    vir_tar_b_batch = [];
                    inf_tar_b_batch = [];

                    if includeAdditionalOutputs ~= 1
                        [isInfected(i)] = POI_BinaryCR_integrated(batchParams(i,:), 0, includeAdditionalOutputs);
                    else
                        [isInfected(i), t_outputs_batch, vir_tar_t_batch, inf_tar_t_batch, vir_tar_b_batch, inf_tar_b_batch] = POI_BinaryCR_integrated(batchParams(i,:), 0, includeAdditionalOutputs);
                    end
    
                    if mod(i,tictocStep) == 0 || i == 1
                        percentageDone = (i/numItersBatch)*100;
                        fprintf('Batch %d of %d, %.2f%% done \n', j, batchNum, percentageDone)
                    end

                    if isInfected(i) == 1 && includeAdditionalOutputs == 1
                        t_outputs{i} = t_outputs_batch';
                        vir_tar_t{i} = vir_tar_t_batch';
                        inf_tar_t{i} = inf_tar_t_batch';
                        vir_tar_b{i} = vir_tar_b_batch';
                        inf_tar_b{i} = inf_tar_b_batch';
                        sim_id{i} = i.*ones(size(t_outputs{i},1),1);
                    end
                end
    
                writeResultsBinaryPOI(batchParams,isInfected,filenameResults,0)

                if includeAdditionalOutputs == 1
                    t_outputs = cell2mat(t_outputs(~cellfun('isempty',t_outputs)));
                    vir_tar_t = cell2mat(vir_tar_t(~cellfun('isempty',vir_tar_t)));
                    inf_tar_t = cell2mat(inf_tar_t(~cellfun('isempty',inf_tar_t)));
                    vir_tar_b = cell2mat(vir_tar_b(~cellfun('isempty',vir_tar_b)));
                    inf_tar_b = cell2mat(inf_tar_b(~cellfun('isempty',inf_tar_b)));
                    sim_id = cell2mat(sim_id(~cellfun('isempty',sim_id)));
    
                    writeResultsAdditionalOutputs(filenameAdditionalOutputs,0,'sim_id',sim_id,'t',t_outputs,'vir_tar_t',vir_tar_t,'inf_tar_t',inf_tar_t,'vir_tar_b',vir_tar_b,'inf_tar_b',inf_tar_b)
                end
    
                batchTime = toc(start)/60;
                cumTime = cumTime + batchTime;
                fprintf(fileID,'Batch %d Simulation Time = %.2f mins\nTotal Running Simulation Time = %.2f mins\n', j, batchTime, cumTime);
                fprintf(fileID,'Number of Infections = %.0f \nNumber of Simulations = %.0f\n\n', sum(isInfected), length(isInfected));

                fprintf('Batch %d Simulation Time = %.2f mins\nTotal Running Simulation Time = %.2f mins\n', j, batchTime, cumTime);
                fprintf('Number of Infections = %.0f \nNumber of Simulations = %.0f\n\n', sum(isInfected), length(isInfected));

                simulationStatusTable = array2table([combination,j]);
                simulationStatusTable.Properties.VariableNames = {'combination','batch'};
                writetable(simulationStatusTable,filenameSimStatus,'WriteMode','overwritesheet','WriteRowNames',true)
            end
    
            if runCluster == 0
                plotparamhist(newparams,'New Parameters')
            end

        end

    fclose(fileID);   

    elseif simConditions == 3

        for combination = 1:size(paramMeanCombinations,1)
            paramsFile = strsplit(filenameParams,'.');
            resultsFile = strsplit(filenameResults,'.');
            paramsFile = strsplit(paramsFile{1},'_');
            resultsFile = strsplit(resultsFile{1},'_');

            completeParams = readtable(sprintf("%s_%d.xlsx", paramsFile{1}, combination));
            if runCluster == 0
                plotparamhist(completeParams,'Complete Parameters')
            end
            resultsComplete = readtable(sprintf("%s_%d.xlsx", resultsFile{1}, combination));
    
            positiveInfection = find(resultsComplete.isInfected==1);
    
            if ~isempty(positiveInfection)
                
                fprintf('Parameters:\n')
                paramsEVAL = completeParams(positiveInfection(1),:)
                makePlots = 1;
                if includeAdditionalOutputs ~= 1
                    [isInfected] = POI_BinaryCR_integrated(paramsEVAL, makePlots, includeAdditionalOutputs);
                else
                    [isInfected, t_outputs, vir_tar_t, inf_tar_t, vir_tar_b, inf_tar_b] = POI_BinaryCR_integrated(paramsEVAL, makePlots, includeAdditionalOutputs);
                end

            else
                fprintf('\n\nNo Infections \n\n')
            end

            errorSim = find(resultsComplete.isInfected==-1);

            if ~isempty(errorSim)
                
                fprintf('Parameters:\n')
                paramsEVAL = completeParams(errorSim(1),:)
                makePlots = 1;
                %[isInfected] = POI_BinaryCR_integrated(paramsEVAL, makePlots);
                if includeAdditionalOutputs ~= 1
                    [isInfected] = POI_BinaryCR_integrated(paramsEVAL, makePlots, includeAdditionalOutputs);
                else
                    [isInfected, t_outputs, vir_tar_t, inf_tar_t, vir_tar_b, inf_tar_b] = POI_BinaryCR_integrated(paramsEVAL, makePlots, includeAdditionalOutputs);
                end
            else
                fprintf('\n\nNo Error Simulations \n\n')
            end

        end

    elseif simConditions == 4

        T = TCell;
        V = Virus;
        I = ICell;
        U = Chemokine;
        D = Drug;
        Geometry = Vagina;
        Mesh = Discretization1D;

        SA = 100; % surface area vaginal epithelium in humans
        Vg = 2; % mL volume of drug gel
        Vs = 2; % mL volume of semen
        Vf = 1; % mL volume of vaginal fluid
        
        % CLOSE TO INFECTION

        % V_0 = 1e5;
        % T0_B = 8e3;
        % T0_T = 3e3;
        % h_E = 0.015;
        % rho = 1400/(24*3600);
        % k_B = 11/(24*3600);
        % k_L = 10/(4*3600);
        % beta = (4e-6)/(24*3600);
        % w = (1e-5)/(24*3600);
        % C0 = 0;
        % T_VD = 0;

        % +VE INFECTION

        % V_0 = 50000.1216804234;
        % h_E = 0.0262852173052513;
        % rho = 0.0328390196877809;
        % k_B = 0.0000362485484221441;
        % k_L = 0.00033704106794302;
        % beta = 2.40847326791511E-11;
        % T0_B = 5837.72001597681;
        % T0_T = 932.455491263257;
        % 
        % w = (1e-5)/(24*3600);
        % C0 = 1e7;
        % T_VD = -1*60*60;

        % EXAMPLE WEAK +VE INFECTION

        % V_0 = 5000;
        % h_E = 0.02;
        % rho = 0.030;
        % k_B = 0.000035;
        % k_L = 0.00035;
        % beta = 2.4E-11;
        % T0_B = 5000;
        % T0_T = 1000;
        % 
        % w = (1e-5)/(24*3600);
        % C0 = 1e7;
        % T_VD = -2*60*60;

        % EXAMPLE STRONG +VE INFECTION TISSUE FIRST

        % V_0 = 4000;
        % h_E = 0.02;
        % rho = 0.025;
        % k_B = 0.00005;
        % k_L = 0.0002;
        % beta = 4E-11;
        % T0_B = 4000;
        % T0_T = 2000;
        % 
        % w = (1e-5)/(24*3600);
        % C0 = 0;
        % T_VD = 0;

        % EXAMPLE STRONG +VE INFECTION BLOOD FIRST
        % V_0 = 5000;
        % h_E = 0.02;
        % rho = 0.010;
        % k_B = 0.00008;
        % k_L = 0.0002;
        % beta = 2E-10;
        % T0_B = 5000;
        % T0_T = 1000;
        % 
        % w = (1e-5)/(24*3600);
        % C0 = 0;
        % T_VD = 0;

        % EXAMPLE STRONG +VE INFECTION BLOOD FIRST v2
        % V_0 = 5000;
        % h_E = 0.02;
        % rho = 0.050;
        % k_B = 0.00008;
        % k_L = 0.0002;
        % beta = 1.5E-11;
        % T0_B = 5000;
        % T0_T = 1000;
        % 
        % w = (1e-5)/(24*3600);
        % C0 = 0;
        % T_VD = 0;

        % EXAMPLE STRONG +VE INFECTION BLOOD FIRST v3
        V_0 = 5000;
        h_E = 0.02;
        rho = 0.050;
        k_B = 0.0001;
        k_L = 0.0002;
        beta = 1E-11;
        T0_B = 10000;
        T0_T = 1000;

        w = (1e-5)/(24*3600);
        C0 = 0;
        T_VD = 0;

        % EXAMPLE STRONG +VE INFECTION BLOOD FIRST v4

        V_0 = 5000;
        h_E = 0.02;
        rho = 0.01;
        k_B = 0.001;
        k_L = 0.0001;
        beta = 3E-11;
        T0_B = 10000;
        T0_T = 1000;

        w = (2e-5)/(24*3600);

        %C0 = 1e7;
        %drugType = 1;
        %T_VD = 0;

        C0 = 30000000;
        drugType = 2;
        T_VD = 0;

        %C0 = 0;
        %drugType = 0;
        %T_VD = 0;

        % EXAMPLE -VE INFECTION

        % V_0 = 5000;
        % h_E = 0.02;
        % rho = 0.010;
        % k_B = 0.0001;
        % k_L = 0.0004;
        % beta = 1E-11;
        % T0_B = 5000;
        % T0_T = 1000;
        % w = (1e-5)/(24*3600);
        % %C0 = 0;
        % %T_VD = 0;
        % C0 = 1e7;
        % T_VD = 0;

        % V_0 = 5000;
        % h_E = 0.02;
        % rho = 0.010;
        % k_B = 0.00008;
        % k_L = 0.0002;
        % beta = 2E-10;
        % T0_B = 5000;
        % T0_T = 1000;
        % 
        % w = (1e-5)/(24*3600);
        % C0 = 0;
        % T_VD = 0;

        % V_0 = 1e4;
        % T0_B = 5e3;
        % T0_T = 5e3;
        % h_E = 0.02;
        % rho = 850/(24*3600);
        % k_B = 3/(24*3600);
        % k_L = 40/(4*3600);
        % beta = (1e-7)/(24*3600);
        % w = (2e-4)/(24*3600);
        % C0 = 1e7;
        % T_VD = -2*60*60;
        %C0 = 0;
        %T_VD = 0;

        paramsEVAL = table(V_0,SA,h_E,rho,k_B,k_L,beta,w,T0_B,T0_T,Vg,Vs,Vf,C0,T_VD,drugType);
        %paramsEVAL = table(V_0,h_E,beta,rho,k_B,T_VD);
        includeAdditionalOutputs = 1;
        makePlots = 1;
        tic
        %[infected,t_inf] = POI_BinaryCR_integrated(paramsEVAL, makePlots)
        
        
        %[isInfected, t_outputs, vir_tar_t, inf_tar_t, vir_tar_b, inf_tar_b] = POI_BinaryCR_integrated(paramsEVAL, makePlots, includeAdditionalOutputs);
        if includeAdditionalOutputs ~= 1
            [isInfected] = POI_BinaryCR_integrated(paramsEVAL, makePlots, includeAdditionalOutputs);
        else
            [isInfected,t_outputs, vir_tar_t, inf_tar_t, vir_tar_b, inf_tar_b] = POI_BinaryCR_integrated(paramsEVAL, makePlots, includeAdditionalOutputs);
        end
        fprintf('isInfected = %d \n', isInfected);
       
        %POI_LumenTissueBlood(paramsEVAL, makePlots)
        
        %infected = POI_BinaryCR_CollTest(paramsEVAL)
        toc
    elseif simConditions == 5

        fileNum = 2;
        paramRowNum = 88;

        T = TCell;
        V = Virus;
        I = ICell;
        U = Chemokine;
        Geometry = Vagina;
        Mesh = Discretization1D;
        
        paramsFile = strsplit(filenameParams,'.');
        resultsFile = strsplit(filenameResults,'.');
        paramsFile = strsplit(paramsFile{1},'_');
        resultsFile = strsplit(resultsFile{1},'_');

        completeParams = readtable(sprintf("%s_%d.xlsx", paramsFile{1}, fileNum));
        paramRow = completeParams(paramRowNum,:)

        makePlots = 1;
        tic
        
        %[isInfected] = POI_BinaryCR_integrated(paramRow, makePlots, includeAdditionalOutputs);
        if includeAdditionalOutputs ~= 1
            [isInfected] = POI_BinaryCR_integrated(paramsEVAL, makePlots, includeAdditionalOutputs);
        else
            [isInfected,t_outputs, vir_tar_t, inf_tar_t, vir_tar_b, inf_tar_b] = POI_BinaryCR_integrated(paramsEVAL, makePlots, includeAdditionalOutputs);
        end
        fprintf('isInfected = %d \n', isInfected);
       
        toc
                
    end
    

