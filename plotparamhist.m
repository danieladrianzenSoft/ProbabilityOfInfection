function plotparamhist(params,title) 
    showFullTitle = 0;
    color_key = 5;
 % Set colors for up to 7 plots
    % color key
    beige    = [149,100,50]/255;
    coral    = [238,0,0]/255;
    grass    = [13,168,42]/255;
    lavender = [139,73,104]/255;
    seafoam  = [67,150,128]/255;
    sky      = [45,131,72]/255;
    steel    = [119,108,94]/255;
    peach    = [253,142,45]/255;
    mustard  = [145,146,39]/255;
    
    map_color = [beige; coral; grass; lavender; seafoam; sky; steel; peach; mustard];

    figure()
    %sgtitle(title, 'FontSize', 28, 'FontWeight','Bold')

    subplot(3,2,1)
    histogram(params.V_0, 'FaceColor', map_color(color_key(1),:), 'FaceAlpha', 0.5, 'EdgeAlpha', 0)
    set(gca, 'FontSize',18)
    %plotHistOverlay('testing','float','off',4,params.V_0)
    if (showFullTitle == 1)
        xlabel('Initial Viral Load, V_{0} (virions/ml)', 'FontSize', 22, 'FontWeight','Bold')
    else
        xlabel('V_{0} (virions/ml)', 'FontSize', 22, 'FontWeight','Bold')
    end

    subplot(3,2,2)
    histogram(params.h_E, 'FaceColor', map_color(color_key(1),:), 'FaceAlpha', 0.5, 'EdgeAlpha', 0)
    set(gca, 'FontSize',18)
    if (showFullTitle == 1)
        xlabel('Epithelial Thickness, h_{E} (cm)', 'FontSize', 22, 'FontWeight','Bold')
    else
        xlabel('h_{E} (cm)', 'FontSize', 22, 'FontWeight','Bold')
    end

    subplot(3,2,3)
    histogram(params.beta*24*3600, 'FaceColor', map_color(color_key(1),:), 'FaceAlpha', 0.5, 'EdgeAlpha', 0)
    set(gca, 'FontSize',18)
    if (showFullTitle == 1)
        xlabel('V infectivity, \beta (cm^3 virion^{-1} day^{-1})', 'FontSize', 22, 'FontWeight','Bold')
    else
        xlabel('\beta (cm^3 virion^{-1} day^{-1})', 'FontSize', 22, 'FontWeight','Bold')
    end

    subplot(3,2,4)
    histogram(params.rho*24*3600, 'FaceColor', map_color(color_key(1),:), 'FaceAlpha', 0.5, 'EdgeAlpha', 0)
    set(gca, 'FontSize',18)
    if (showFullTitle == 1)
        xlabel('Production Rate of V from I, \rho (day^{-1})', 'FontSize', 22, 'FontWeight','Bold')
    else
        xlabel('\rho (day^{-1})', 'FontSize', 22, 'FontWeight','Bold')
    end

    subplot(3,2,5)
    histogram(params.k_B*24*3600, 'FaceColor', map_color(color_key(1),:), 'FaceAlpha', 0.5, 'EdgeAlpha', 0)
    set(gca, 'FontSize',18)
    if (showFullTitle == 1)
        xlabel('Clearance Rate of V from Stroma to Blood, k_{B} (day^{-1})', 'FontSize', 22, 'FontWeight','Bold')
    else
        xlabel('k_{B} (day^{-1})', 'FontSize', 22, 'FontWeight','Bold')
    end

    subplot(3,2,6)
    histogram(params.k_L*24*3600, 'FaceColor', map_color(color_key(1),:), 'FaceAlpha', 0.5, 'EdgeAlpha', 0)
    set(gca, 'FontSize',18)
    if (showFullTitle == 1)
        xlabel('Clearance Rate of V from Blood, c (day^{-1})', 'FontSize', 22, 'FontWeight','Bold')
    else
        xlabel('c (day^{-1})', 'FontSize', 22, 'FontWeight','Bold')
    end






end