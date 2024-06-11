% Specify the folder where the files are located
folderPath = './Cluster/Cluster_03122024_CellCellVirCellOutput_1';

% List all Excel files starting with 'AdditionalOutputs'
files = dir(fullfile(folderPath, 'AdditionalOutputs*.xlsx'));

% List column names in excel file
colNames = {'vir_tar_t', 'inf_tar_t', 'vir_tar_b', 'inf_tar_b'}; % Column names

% Loop through each file
for k = 1:length(files)
    % Read the current file
    fileName = fullfile(files(k).folder, files(k).name);
    opts = detectImportOptions(fileName);
    data = readtable(fileName, opts);

    % Extract unique time points from this file
    timePoints = unique(data.t);

    % Initialize containers to store full datasets for each variable across all files
    fullDataSets = struct('vir_tar_t', [], 'inf_tar_t', [], 'vir_tar_b', [], 'inf_tar_b', []);

    % Accumulate data for full datasets
    for j = 1:numel(colNames)
        if isfield(fullDataSets, colNames{j})
            fullDataSets.(colNames{j}) = [fullDataSets.(colNames{j}), data{:, colNames{j}}];
        else
            fullDataSets.(colNames{j}) = data{:, colNames{j}};
        end
    end

    % Initialize matrices to hold mean (or median) values and confidence intervals for this file
    means = zeros(length(timePoints), numel(colNames)); % Change 'means' to 'medians' if using median
    medians = zeros(length(timePoints), numel(colNames)); % Change 'means' to 'medians' if using median

    confIntervals = zeros(length(timePoints), numel(colNames), 2); % Lower and upper bounds

    % Calculate mean/median and confidence intervals for each time point
    for i = 1:length(timePoints)
        t = timePoints(i);
        for j = 1:numel(colNames) % Loop through the four columns of interest
            colName = colNames{j}; % Access column name
            dataAtTimePoint = data{data.t == t, colName};

            % Calculate mean or median
            means(i, j) = mean(dataAtTimePoint); % Replace with median(dataAtTimePoint) if you prefer median
            medians(i, j) = median(dataAtTimePoint);

            % Calculate confidence interval (assuming normal distribution)
            SEM = std(dataAtTimePoint) / sqrt(length(dataAtTimePoint)); % Standard Error
            ts = tinv([0.025  0.975], length(dataAtTimePoint)-1); % T-Score for 95% CI
            confIntervals(i, j, :) = means(i, j) + ts*SEM; % Adjust for median if necessary
        end
    end

    % Here, 'means' contains the mean values for each column at each time point for the current file,
    % and 'confIntervals' contains the corresponding 95% confidence intervals.
    % You can now save or plot these results as needed, file by file.
    % For example, to save the results to a .mat file:
    %save(fullfile(folderPath, sprintf('Results_%s.mat', files(k).name)), 'means', 'confIntervals', 'timePoints');

    % Convert fullDataSets struct to cell array format expected by createAdditionalOutputPlots
    fullDataSetsArray = {fullDataSets.vir_tar_t, fullDataSets.inf_tar_t, fullDataSets.vir_tar_b, fullDataSets.inf_tar_b};

    % Or, to simply display which file has been processed:
    fprintf('Processed file: %s\n', files(k).name);
    

    additionalOutputSettingLabels = ["showCellCellAndVirusCellAllCompartments", "showCellCellAndVirusCellPerCompartment", "showFullDataSets"];
    if k == 1
        additionalOutputSettingValues = [1,1,1];
    else
        additionalOutputSettingValues = [0,0,0];
    end

    additionalOutputSettings = dictionary(additionalOutputSettingLabels, additionalOutputSettingValues);
    createAdditionalOutputPlots(timePoints/3600,means(:,1),means(:,2),means(:,3),means(:,4),additionalOutputSettings,confIntervals,fullDataSetsArray)
end

% % Specify the folder where the files are located
% folderPath = './Cluster/Cluster_03122024_CellCellVirCellOutput_1';
% 
% % List all Excel files starting with 'AdditionalOutputs'
% files = dir(fullfile(folderPath, 'AdditionalOutputs*.xlsx'));
% 
% % Variables to process
% colNames = {'vir_tar_t', 'inf_tar_t', 'vir_tar_b', 'inf_tar_b'}; % Column names
% 
% % Loop through each file
% for k = 1:length(files)
%     % Read the current file
%     fileName = fullfile(files(k).folder, files(k).name);
%     opts = detectImportOptions(fileName);
%     data = readtable(fileName, opts);
% 
%     % Extract unique time points from this file
%     timePoints = unique(data.t);
% 
%     % Initialize matrices to hold mean values and confidence intervals for this file
%     means = zeros(length(timePoints), numel(colNames)); % Mean values
%     confIntervals = zeros(length(timePoints), numel(colNames), 2); % Confidence intervals
% 
%     % Plot settings
%     lightGray = [0.8, 0.8, 0.8]; % Light gray for individual series
%     black = [0, 0, 0]; % Black for mean
%     red = [1, 0, 0]; % Red for confidence intervals
% 
%     % For each variable, process data and plot
%     for j = 1:numel(colNames)
%         varName = colNames{j};
%         figure; hold on; % Create a new figure for each variable
%         title([varName ' across all simulations'], 'Interpreter', 'none');
%         xlabel('Time');
%         ylabel(varName);
% 
%         % Plot all individual time series in light gray
%         for m = 1:size(data, 1)
%             plot(data.t(m), data{m, varName}, 'Color', lightGray, 'LineWidth', 0.5);
%         end
% 
%         % Calculate mean, confidence intervals, and plot them
%         for i = 1:length(timePoints)
%             t = timePoints(i);
%             dataAtTimePoint = data{data.t == t, varName};
% 
%             % Calculate mean
%             means(i, j) = mean(dataAtTimePoint);
% 
%             % Calculate confidence interval (assuming normal distribution)
%             SEM = std(dataAtTimePoint) / sqrt(length(dataAtTimePoint)); % Standard Error
%             ts = tinv([0.025  0.975], length(dataAtTimePoint)-1); % T-Score for 95% CI
%             confIntervals(i, j, :) = means(i, j) + ts*SEM; % Adjust for median if necessary
%         end
% 
%         % Plot mean in black
%         plot(timePoints, means(:, j), 'Color', black, 'LineWidth', 2);
% 
%         % Plot confidence intervals in red
%         ciUpper = squeeze(confIntervals(:, j, 2));
%         ciLower = squeeze(confIntervals(:, j, 1));
%         fill([timePoints; flipud(timePoints)], [ciLower; flipud(ciUpper)], red, 'FaceAlpha', 0.1, 'EdgeColor', 'none');
% 
%         legend({'Individual Simulations', 'Mean', '95% Confidence Interval'}, 'Location', 'Best');
%         hold off;
%     end
% 
%     % Indicate completion of processing for this file
%     fprintf('Processed and plotted file: %s\n', files(k).name);
% end