clear
clc

B_POI = readtable('B_POI.xlsx');

numInfected = sum(B_POI.infected);

POI = numInfected/size(B_POI,1)*100

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

% noDrug = readtable('B_POI_0916.xlsx');
% drug0hrs = readtable('B_POI_0917a.xlsx');
% drugn2hrs = readtable('B_POI_0917b.xlsx');
% drugp2hrs = readtable('B_POI_0918.xlsx');
% 
% noDrugPOI = sum(noDrug.infected)/length(noDrug.infected);
% drug0hrsPOI = sum(drug0hrs.infected)/length(drug0hrs.infected);
% drugn2hrsPOI = sum(drugn2hrs.infected)/length(drugn2hrs.infected);
% drugp2hrsPOI = sum(drugp2hrs.infected)/length(drugp2hrs.infected);

%POI = [0.255,0.16,0.15,0.16];
%POI = [noDrugPOI, drug0hrsPOI, drugn2hrsPOI, drugp2hrsPOI]*100;
labels = categorical({'No Drug'});
%labels = reordercats(labels,{'No Drug','Drug, 0rs delay','Drug, 2hrs prior','Drug, 2hrs post'});

% Y = [10 21 33 52];
% bar(X,Y)
bar(labels,POI,'FaceColor', map_color(color_key(1),:),'FaceAlpha', 0.5, 'EdgeAlpha', 0);
set(gca, 'FontSize', 18)
ylabel('Probability of Infection (%)', 'FontSize',22,'FontWeight','Bold')

%B_POI.t_inf(noDrug.t_inf>0)

infectionTime = B_POI.t_inf(B_POI.t_inf>0)/(24*60*60);

histogram(infectionTime, 'FaceColor', map_color(color_key(1),:), 'FaceAlpha', 0.5, 'EdgeAlpha', 0)
xlabel('Time of first infected cell (days)', 'FontSize',16,'FontWeight','Bold')
ylabel('Frequency', 'FontSize',16,'FontWeight','Bold')
