function createAdditionalOutputPlots(t, vir_tar_t, inf_tar_t, vir_tar_b, inf_tar_b, settings, varargin)
    lineWidth = 6;
    if t(end)/24 > 18
        tend = 18;
    else
        tend = t(end)/24;
    end

    % Check if confidence intervals are provided
    if ~isempty(varargin) && length(varargin) == 1
        confIntervals = varargin{1}; % Expecting {fullDataSet_vir_tar_t, fullDataSet_inf_tar_t, fullDataSet_vir_tar_b, fullDataSet_inf_tar_b}
        fullDataSets = {};
    end
    % Check if full datasets are provided
    if ~isempty(varargin) && length(varargin) > 1
        confIntervals = varargin{1};
        fullDataSets = varargin{2}; % Expecting {fullDataSet_vir_tar_t, fullDataSet_inf_tar_t, fullDataSet_vir_tar_b, fullDataSet_inf_tar_b}
    end
    %Check if confIntervals are provided
    if isempty(varargin)
        confIntervals = {};
        fullDataSets = {};
    end

    if isKey(settings,"showCellCellAndVirusCellAllCompartments") && settings("showCellCellAndVirusCellAllCompartments") == 1
        f = figure();
        f.Position = [383, 363, 793, 627];
        ax = axes('Parent', f);
        cm = getCustomColormap(6); % Assuming getCustomColormap is defined elsewhere
        
        % Plot the data
        plot(ax, t/24, vir_tar_t, 'color', cm(1,:), 'LineWidth', lineWidth);
        hold on;
        plot(ax, t/24, inf_tar_t, 'color', cm(2,:), 'LineWidth', lineWidth);
        plot(ax, t/24, vir_tar_b, 'color', cm(3,:), 'LineWidth', lineWidth);
        plot(ax, t/24, inf_tar_b, 'color', cm(4,:), 'LineWidth', lineWidth);

        % Plot confidence intervals if provided
        if ~isempty(confIntervals)
            %confIntervals = varargin{1}; % 3D matrix of confidence intervals
            % Adjust indices based on whether we're plotting tissue or blood data

            % Plot CI for data1
            ci_data1 = squeeze(confIntervals(:, 1, :));
            fill([t/24; flipud(t/24)], [ci_data1(:,1); flipud(ci_data1(:,2))], cm(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

            % Plot CI for data2
            ci_data2 = squeeze(confIntervals(:, 2, :));
            fill([t/24; flipud(t/24)], [ci_data2(:,1); flipud(ci_data2(:,2))], cm(2,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

            % Plot CI for data3
            ci_data3 = squeeze(confIntervals(:, 3, :));
            fill([t/24; flipud(t/24)], [ci_data3(:,1); flipud(ci_data3(:,2))], cm(3,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');

            % Plot CI for data4
            ci_data4 = squeeze(confIntervals(:, 4, :));
            fill([t/24; flipud(t/24)], [ci_data4(:,1); flipud(ci_data4(:,2))], cm(4,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
        end

        % Customize the plot
        xlabel(ax, 'Time (days)', 'FontSize', 40, 'FontWeight', 'Bold');
        ylabel(ax, '# cell infections / mL', 'FontSize', 40, 'FontWeight', 'Bold');
        title(ax, 'Cell infections by source and location', 'FontSize', 44, 'FontWeight', 'Bold');
        legend(ax, {'Virus-cell in tissue', 'Cell-cell in tissue','Virus-cell in blood', 'Cell-cell in blood'}, "Location", "northwest");
        legend boxoff;
        set(ax, 'FontSize', 40);
        hold off;
    

    end
    if isKey(settings,"showCellCellAndVirusCellPerCompartment") && settings("showCellCellAndVirusCellPerCompartment") == 1
        % Define the names and labels for easier reference and modification
        dataPairs = {...
        {vir_tar_t, inf_tar_t, 'Free virus-cell', 'Cell-cell', 'Cell Infections in Tissue', '# infections / mL'}, ...
        {vir_tar_b, inf_tar_b, 'Free virus-cell', 'Cell-cell', 'Cell Infections in Blood', '# infections / mL'}};

    
        % Iterate through each data pair for plotting
        for i = 1:length(dataPairs)
            data1 = dataPairs{i}{1};
            data2 = dataPairs{i}{2};
            data1Name = dataPairs{i}{3};
            data2Name = dataPairs{i}{4};
            titleStr = dataPairs{i}{5};
            yLabel = dataPairs{i}{6};
            
            % Create figure and axes
            f = figure();
            f.Position = [383, 363, 793, 627];
            ax = axes('Parent', f);
            cm = getCustomColormap(3); % Assuming getCustomColormap is defined elsewhere
    
            % Plot the data
            plot(ax, t/24, data1, 'color', cm(1,:), 'LineWidth', lineWidth);
            hold on;
            plot(ax, t/24, data2, 'color', cm(2,:), 'LineWidth', lineWidth);
    
            % Plot confidence intervals if provided
            if ~isempty(confIntervals)
                %confIntervals = varargin{1}; % 3D matrix of confidence intervals
                % Adjust indices based on whether we're plotting tissue or blood data
                indexOffset = (i-1)*2;

                % Plot CI for data1
                ci_data1 = squeeze(confIntervals(:, indexOffset+1, :));
                fill([t/24; flipud(t/24)], [ci_data1(:,1); flipud(ci_data1(:,2))], cm(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
    
                % Plot CI for data2
                ci_data2 = squeeze(confIntervals(:, indexOffset+2, :));
                fill([t/24; flipud(t/24)], [ci_data2(:,1); flipud(ci_data2(:,2))], cm(2,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
            end

            % Customize the plot
            xlabel(ax, 'Time (days)', 'FontSize', 40, 'FontWeight', 'Bold');
            ylabel(ax, yLabel, 'FontSize', 40, 'FontWeight', 'Bold');
            title(ax, titleStr, 'FontSize', 44, 'FontWeight', 'Bold');
            legend(ax, {data1Name, data2Name}, "Location", "northwest");
            legend boxoff;
            set(ax, 'FontSize', 40);
            hold off;
        end

        if isKey(settings,"showFullDataSets") && settings("showFullDataSets") == 1
            if ~isempty(fullDataSets)
                numTimesteps = length(t);
                for dsIndex = 1:length(fullDataSets)
                    switch dsIndex
                        case 1
                            meanToPlot = vir_tar_t;
                            ciToPlot = squeeze(confIntervals(:, 1, :));
                            yLabel = "# infections / mL";
                            titleStr = "Virus-cell infections in tissue";
                        case 2
                            meanToPlot = inf_tar_t;
                            ciToPlot = squeeze(confIntervals(:, 2, :));
                            yLabel = "# infections / mL";
                            titleStr = "Cell-cell infections in tissue";
                        case 3
                            meanToPlot = vir_tar_b;
                            ciToPlot = squeeze(confIntervals(:, 3, :));
                            yLabel = "# infections / mL";
                            titleStr = "Virus-cell infections in blood";
                        case 4
                            meanToPlot = inf_tar_b;
                            ciToPlot = squeeze(confIntervals(:, 4, :));
                            yLabel = "# infections / mL";
                            titleStr = "Cell-cell infections in blood";
                    end
                    f = figure();
                    f.Position = [383, 363, 793, 627];
                    ax = axes('Parent', f);
                    currentOutput = fullDataSets{dsIndex};
                    numSeries = length(fullDataSets{dsIndex})/numTimesteps;
                    for seriesIndex = 1:numSeries
                        range = (seriesIndex-1)*numTimesteps+1:seriesIndex*numTimesteps;
                        plot(ax, t/24, currentOutput(range), 'Color', [0.8, 0.8, 0.8], 'LineWidth', 1);
                        hold on;
                    end
    
                    plot(ax, t/24, meanToPlot, 'Color', [0, 0, 0], 'LineWidth', 3);
                    fill([t/24; flipud(t/24)], [ciToPlot(:,1); flipud(ciToPlot(:,2))], cm(1,:), 'FaceAlpha', 0.3, 'EdgeColor', 'none');
                    
                    xlabel(ax, 'Time (days)', 'FontSize', 40, 'FontWeight', 'Bold');
                    ylabel(ax, yLabel, 'FontSize', 40, 'FontWeight', 'Bold');
                    title(ax, titleStr, 'FontSize', 44, 'FontWeight', 'Bold');
                    %legend(ax, {data1Name, data2Name}, "Location", "northwest");
                    %legend boxoff;
                    set(ax, 'FontSize', 40);
                    hold off;
    
                end
            end
        end

    end

end