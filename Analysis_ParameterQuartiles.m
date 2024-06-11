clear
clc

parameterOfInterest = 'k_B';

% Specify the folder where the files are located
folderPath = './';

% List all Excel files starting with 'resultsPOI'
files = dir(fullfile(folderPath, 'resultsPOI*.xlsx'));

% Initialize an empty table for aggregation with predefined columns if known
% If the structure of all files is consistent, you can initialize the table with column names
%columnsLabels = {'V_0', 'SA', 'h_E', 'rho', 'k_B', 'k_L', 'beta', 'w', 'T0_B', 'T0_T', 'Vf', 'Vs', 'Vg', 'C_G0', 'T_VD'};
%allData = table('Size', [0, length(columnsLabels)], 'VariableNames', columnLabels);

% Initialize a cell array to hold data from each file
dataCells = cell(1, length(files));

for k = 1:length(files)
    % Read the current file
    fileName = fullfile(files(k).folder, files(k).name);
    opts = detectImportOptions(fileName);
    data = readtable(fileName, opts);
    
    % Store the table in the cell array
    dataCells{k} = data;
end

% Concatenate all tables from the cell array into one table
allData = vertcat(dataCells{:});

% Convert parameter of interest to numeric only if it is not already numeric
if iscell(allData.(parameterOfInterest))
    allData.(parameterOfInterest) = str2double(allData.(parameterOfInterest));
end

% Calculate the quartiles for parameter of interest
quartiles = quantile(allData.(parameterOfInterest), [0.25 0.5 0.75]);

% Add a new column to classify each row into a quartile group
allData.quartileGroup = zeros(height(allData), 1); % Initialize column
for i = 1:height(allData)
    if allData.(parameterOfInterest)(i) <= quartiles(1)
        allData.quartileGroup(i) = 1; % First quartile
    elseif allData.(parameterOfInterest)(i) <= quartiles(2)
        allData.quartileGroup(i) = 2; % Second quartile
    elseif allData.(parameterOfInterest)(i) <= quartiles(3)
        allData.quartileGroup(i) = 3; % Third quartile
    else
        allData.quartileGroup(i) = 4; % Fourth quartile
    end
end

% Now, 'allData' contains an additional column 'quartileGroup' 
% that indicates the quartile group for each row based on parameter of
% interest

num_phases = 3;
num_quartiles = 4;
num_batches = 5;

n = zeros(num_quartiles,num_phases);
batch_sizes = zeros(num_quartiles*num_batches,num_phases);
n_infections = zeros(num_quartiles,num_phases);
n_infections_batched = zeros(num_quartiles*num_batches,num_phases);

v0_mean = zeros(num_quartiles,1);
v0_median = zeros(num_quartiles,1);
v0_variance = zeros(num_quartiles,1);

for q = 1:num_quartiles
    % Extract data for the current quartile group
    quartileData = allData(allData.quartileGroup == q, :);

    for p = 1:num_phases
        %PHASE = 0: MIDCYCLE, 1: FOLLICULAR, 2: LUTEAL
        quartileDataPhase = quartileData(quartileData.phase == (p-1), :);
        n(q,p) = height(quartileDataPhase);
        n_infections(q,p) = sum(quartileDataPhase.isInfected);

        batchSize = floor(n(q,p)/num_batches);

        for b = 1:num_batches
            if b < num_batches
               batch_indices = batchSize*(b-1) + 1 : batchSize*b;
            else
               batch_indices =  batchSize*(b-1) + 1 : n(q,p);
            end
            quartileDataPhase_batched = quartileDataPhase(batch_indices,:);
            batchIndex = (q - 1) * num_batches + b; % This creates a unique index for every batch across all quartiles and phases

            batch_sizes(batchIndex,p) = length(batch_indices);
            n_infections_batched(batchIndex,p) = sum(quartileDataPhase_batched.isInfected);
        end
    end

    v0_mean(q) = mean(quartileData.(parameterOfInterest));
    v0_median(q) = median(quartileData.(parameterOfInterest));
    v0_variance(q) = var(quartileData.(parameterOfInterest));

end

poi = n_infections/n;
poi_batched = n_infections_batched ./ batch_sizes;

% Reshape the matrix into a 3D matrix where the third dimension is quartiles
poi_batched_3D = reshape(poi_batched', [num_phases, num_batches, num_quartiles]);
batch_sizes_3D = reshape(batch_sizes', [num_phases, num_batches, num_quartiles]);

% Permute the dimensions to [batches, phases, quartiles]
poi_batched_permuted = permute(poi_batched_3D, [2, 1, 3]);
batch_sizes_permuted = permute(batch_sizes_3D, [2, 1, 3]);

% Finally, reshape this 3D matrix into a single column vector
poi_batched_column = reshape(poi_batched_permuted, [], 1);
batch_sizes_column = reshape(batch_sizes_permuted, [], 1);

% Get batch numbers
batches_column = repmat(linspace(1,num_batches,num_batches),[1,num_quartiles*num_phases])';

% Get phase numbers
phases_column = repmat([zeros(num_batches, 1); ones(num_batches, 1); 2*ones(num_batches, 1)], num_quartiles, 1);

% Get quartile numbers
quartiles_num_column = [1.*ones(num_batches*num_phases, 1); 
                        2.*ones(num_batches*num_phases, 1); 
                        3.*ones(num_batches*num_phases, 1); 
                        4.*ones(num_batches*num_phases, 1)];

quartiles_column = 'Q' + string(quartiles_num_column);

% Construct table
column_names = ["Batch", "Batch Size", sprintf('%s Quartile',parameterOfInterest), "Phase", "POI"];
output = table(batches_column, batch_sizes_column, quartiles_column, phases_column, poi_batched_column*100, 'VariableNames', column_names);
