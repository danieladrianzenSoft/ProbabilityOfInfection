    clear
    clc
    
    % simConditions:
    % 0 - overwrite entire dataset.
    % 1 - restart simulation from last saved batch. 
    % 2 - add new parameters to existing simulation, evaluate batches and append.
    % 3 - plot params and evaluate representative result for prior sim where isInfected = 1.
    % 4 - run independent sim with specified parameters.

    simConditions = 4;
    runCluster = 0;

%     N = 20000;
%     tictocStep = 200;
%     batchSize = 1000;
    N = 100;
    tictocStep = 20;
    batchSize=25;
    C0 = 0;
    DVdelay = 0;
    

    % label output files with date stamp
    %dateStamp = sprintf('%s', datestr(now, 'yy/mm/dd-HH:MM'));

    dateStamp = string(datetime("now"));

    filenameParams = 'parameters.xlsx';
    filenameResults = 'B_POI.xlsx';
    filenameConfig = 'SimulationConfig.txt';
    
    if simConditions < 3

        if N > batchSize
            batchNum = floor(N/batchSize);
        else
            batchNum = 1;
        end

        if simConditions == 0
            overwrite = 1;
            newparams = samplingNormalPOI_BinaryCR(filenameParams,N,C0,DVdelay,overwrite);
            OutputPlaceholderTable = array2table(zeros(0,size(newparams,2)+1));
            %OutputPlaceholderTable.Properties.VariableNames = {'V_0','h_E', ...
            %   'rho','k_B','k_L','C_G0','T_VD','infected','t_inf'};
            OutputPlaceholderTable.Properties.VariableNames = {'V_0','h_E', ...
               'rho','k_B','k_L','beta','C_G0','T_VD','isInfected'};
            writetable(OutputPlaceholderTable,filenameResults,'WriteMode','overwritesheet','WriteRowNames',true)
            fileID = fopen(filenameConfig,'w+');
            fprintf(fileID,'\nDate Started: %s\n',dateStamp);
            fprintf(fileID,'Number of Iterations: %d\n',N);
            fprintf(fileID,'Number of Batches: %d\n',batchNum);
            startBatch = 1;
        elseif simConditions == 1
            overwrite = 0;
            newparams = readtable(filenameParams);
            fileID = fopen(filenameConfig,'a');
            startBatch = (size(readtable(filenameResults),1)/batchSize)+1;
            fprintf(fileID,'\n ------------------------------- \n');
            fprintf(fileID,'\nDate Resumed: %s\n',dateStamp);
            fprintf(fileID,'Starting Batch: %d\n',startBatch);
        elseif simConditions == 2
            overwrite = 0;
            newparams = samplingNormalPOI_BinaryCR(filenameParams,N,C0,DVdelay,overwrite);
            completeParams = readtable(filenameParams);
            fileID = fopen(filenameConfig,'a');
            fprintf(fileID,'\n ------------------------------- \n');
            fprintf(fileID,'\nDate Restarted with New Batches: %s\n',dateStamp);
            fprintf(fileID,'Number of Additional Iterations: %d\n',N);
            fprintf(fileID,'Number of Additional Batches: %d\n',batchNum);
            fprintf(fileID,'Total Num ber of Iterations: %d\n',size(completeParams,1));
            startBatch = 1;        
        end

        cumTime = 0;
        
        %plotparamhist(newparams,'New Parameters')
        
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

            %p4 = Par(numIters);

            parfor i = 1:numItersBatch
                [isInfected(i)] = POI_BinaryCR_integrated(batchParams(i,:), 0);

                if mod(i,tictocStep) == 0 || i == 1
                    percentageDone = (i/numItersBatch)*100;
                    fprintf('Batch %d of %d, %.2f%% done \n', j, batchNum, percentageDone)
                end
                %    p4(i) = Par.toc;
            end
            %stop(p4);

            writeResultsBinaryPOI(batchParams,isInfected,filenameResults,0)

            batchTime = toc(start)/60;
            cumTime = cumTime + batchTime;
            fprintf(fileID,'Batch %d Simulation Time = %.2f mins\nTotal Running Simulation Time = %.2f mins\n\n', j, batchTime, cumTime);
            %if runCluster == 0
            fprintf('Batch %d Simulation Time = %.2f mins\nTotal Running Simulation Time = %.2f mins\n\n', j, batchTime, cumTime);
            %end
        end

        if runCluster == 0
            plotparamhist(newparams,'New Parameters')
        end

    fclose(fileID);   

    elseif simConditions == 3  
        
        completeParams = readtable(filenameParams);
        if runCluster == 0
            plotparamhist(completeParams,'Complete Parameters')
        end

        resultsComplete = readtable(filenameResults);

        positiveInfection = find(resultsComplete.isInfected==1);

        if ~isempty(positiveInfection)
            
            fprintf('Parameters:\n')
            paramsEVAL = completeParams(positiveInfection(1),:)
            makePlots = 1;
            [isInfected] = POI_BinaryCR_integrated(paramsEVAL, makePlots);
        else
            fprintf('\n\nNo Infections \n\n')
        end

    elseif simConditions == 4

        T = TCell;
        V = Virus;
        I = ICell;
        U = Chemokine;
        Geometry = Vagina;
        Mesh = Discretization1D;
        
        V_0 = 1e3;
        %V_0 = 1e2;
        h_E = 0.08;
        rho = 850/(24*3600);
        %rho = 1.5e3/(24*3600);
        %rho = 85/(24*3600);
        %rho = 1.5e3/(24*3600);
        %rho = 10/(24*3600);
        %k_B = 3/(24*3600);
        k_B = 5/(24*3600);
        %k_B = 10/(24*3600);

        k_L = 10/(24*3600);
        %C_G0 = 10^7;
        %C_G0 = 0;
        beta = (6.5*10^(-7))/(24*3600);
        %beta = (6.5*10^(-8))/(24*3600);
        %beta = 1e-5/(24*3600);
        %T_VD = 0; 
        paramsEVAL = table(V_0,h_E,rho,k_B,k_L,beta);
        %paramsEVAL = table(V_0,h_E,beta,rho,k_B,T_VD);
        makePlots = 1;
        tic
        %[infected,t_inf] = POI_BinaryCR_integrated(paramsEVAL, makePlots)
        
        [isInfected] = POI_BinaryCR_integrated(paramsEVAL, makePlots);

        fprintf('isInfected = %d \n', isInfected);
       
        %POI_LumenTissueBlood(paramsEVAL, makePlots)
        
        %infected = POI_BinaryCR_CollTest(paramsEVAL)
        toc
        
    end
