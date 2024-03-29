function createHeatmap(t,Mesh,DrugTissue,DrugDPTissue,HIVTissue,drugApplied)

    x = Mesh.x;
    endT = 18; %days
    [XX,YY] = meshgrid(t/24,x);
    %ZZ_V = real(log10(HIVTissue));
    %ZZ_V(ZZ_V<0)=0;
    ZZ_V = HIVTissue;
    %ZZ_V(ZZ_V<0)=0;
    CTopHIV = floor(max(max(ZZ_V))) + floor( max(max(ZZ_V)-floor(max(max(ZZ_V))))/0.5) * 0.5;

    IC50 = 180;
    calcQ = @(c) 1./(1+(IC50./(c)));

    cm = getCustomColormap(800,"colormap","purple-teal");

    indLE = Mesh.numL;
    indES = Mesh.numL+Mesh.numE;
    lineLE_X = [min(t),max(t)];
    lineLE_Y = x(indLE).*ones(2,1);
    lineES_X = [min(t),max(t)];
    lineES_Y = x(indES).*ones(2,1);
    lineColor = [0.84,0.84,0.84];
    xTickVec = [0,5,10,15,20];

      if drugApplied == 1
         f1 = figure();
         f1.Position = [383,363,793,627];
         ax1 = axes('Parent',f1);
         %ZZ_D = log10(DrugTissue);  
         %ZZ_D(ZZ_D<0) = 0;
         ZZ_D = DrugTissue;
         ZZ_D(ZZ_D<0) = 0;
         %ZZ_1 = DrugTissue;
         CTopDrug = floor(max(max(ZZ_D))) + floor( max(max(ZZ_D)-floor(max(max(ZZ_D))))/0.5) * 0.5;
         %ax1 = nexttile;
         %contourf(XX,YY,ZZ_1,'LineColor','none')
         %surf(ax1,XX,YY,ZZ_D,'EdgeColor','none')
         p = pcolor(ax1,XX,YY,ZZ_D);
         set(p, 'EdgeColor', 'none');
         title('TFV Concentration (ng/ml)','FontSize',22,'FontWeight','Bold')
         view([0,90])
         set(gca,'Ydir','reverse')
         xlim([0,endT])
         ylim([0,max(x)])
         clim(ax1,[1e1,1e7])
         xticks(ax1,xTickVec)
         yticks(ax1,[0,0.1,0.2,0.3])
         colormap(ax1,'hot')
         set(gca,'ColorScale','log')
         colorbar('Ticks',[1e1,1e3,1e5,1e7])
         %clim(ax1,[0,CTopDrug])
         xlabel('Time (days)','FontSize',40,'FontWeight','Bold')
         ylabel('Depth (cm)','FontSize',40,'FontWeight','Bold')
         XT = get(gca,'XTick');
         set(gca,'FontSize',40)
         YT = get(gca,'YTick');
         set(gca,'FontSize',40)
         hold on
         plot(ax1,lineLE_X,lineLE_Y,'LineWidth',2,'Color',lineColor,'LineStyle','--')
         plot(ax1,lineES_X,lineES_Y,'LineWidth',2,'Color',lineColor,'LineStyle','--')

         f2 = figure();
         f2.Position = [383,363,793,627];
         ax2 = axes('Parent',f2);
         ZZ_DP = DrugDPTissue;  
         ZZ_DP(ZZ_DP<0) = 0;
         %ZZ_DP = DrugDPTissue;
         %ZZ_DP(ZZ_DP<0) = 0;
         %surf(ax2,XX,YY,ZZ_DP,'EdgeColor','none')
         p = pcolor(ax2,XX,YY,ZZ_DP);
         set(p, 'EdgeColor', 'none');
         title('TFV-DP Concentration (ng/ml)','FontSize',22,'FontWeight','Bold')
         view([0,90])
         set(gca,'Ydir','reverse')
         xlim([0,endT])
         ylim([0,max(x)])
         clim(ax2,[1e-1,1e5])
         xticks(ax2,xTickVec)
         yticks(ax2,[0,0.1,0.2,0.3])
         colormap(ax2,'winter')
         set(gca,'ColorScale','log')
         colorbar(ax2,'Ticks',[1e0,1e2,1e4])
         xlabel('Time (days)','FontSize',40,'FontWeight','Bold')
         ylabel('Depth (cm)','FontSize',40,'FontWeight','Bold')
         XT = get(gca,'XTick');
         set(gca,'FontSize',40)
         YT = get(gca,'YTick');
         set(gca,'FontSize',40)
         hold on
         plot(ax2,lineLE_X,lineLE_Y,'LineWidth',2,'Color',lineColor,'LineStyle','--')
         plot(ax2,lineES_X,lineES_Y,'LineWidth',2,'Color',lineColor,'LineStyle','--')

         f3 = figure();
         f3.Position = [383,363,793,627];
         ax3 = axes('Parent',f3);
         ZZ_QDP = calcQ(DrugDPTissue);  
         %ZZ_QDP(ZZ_QDP<0) = 0;
         %surf(ax3,XX,YY,ZZ_QDP,'EdgeColor','none')
         p = pcolor(ax3,XX,YY,ZZ_QDP);
         set(p, 'EdgeColor', 'none');
         title('Drug efficacy','FontSize',22,'FontWeight','Bold')
         view([0,90])
         set(gca,'Ydir','reverse')
         xlim([0,endT])
         ylim([0,max(x)])
         clim(ax3,[0,1])
         xticks(ax3,xTickVec)
         yticks(ax3,[0,0.1,0.2,0.3])
         colormap(ax3,'parula')
         colorbar(ax3,'Ticks',[0,0.5,1])
         xlabel('Time (days)','FontSize',40,'FontWeight','Bold')
         ylabel('Depth (cm)','FontSize',40,'FontWeight','Bold')
         XT = get(gca,'XTick');
         set(gca,'FontSize',40)
         YT = get(gca,'YTick');
         set(gca,'FontSize',40)
         hold on
         plot(ax3,lineLE_X,lineLE_Y,'LineWidth',2,'Color',lineColor,'LineStyle','--')
         plot(ax3,lineES_X,lineES_Y,'LineWidth',2,'Color',lineColor,'LineStyle','--')

      end

         %ax2 = nexttile;
         f = figure();
         f.Position = [383,363,793,627];
         ax = axes('Parent',f);
         %contourf(XX,YY,ZZ_2,'LineColor','none')
         %surf(ax,XX,YY,ZZ_V,'EdgeColor','none')
         p = pcolor(ax,XX,YY,ZZ_V);
         set(p, 'EdgeColor', 'none');
         title(ax,'Viral Load (virions/ml)','FontSize',22,'FontWeight','Bold')
         view([0,90])
         set(gca,'Ydir','reverse')
         xlim([0,endT])
         ylim([0,max(x)])
         clim(ax,[1e-2,1e5])
         xticks(ax,xTickVec)
         yticks(ax,[0,0.1,0.2,0.3])
         colormap(ax,'pink')
         colorbar(ax,'Ticks',[1e-2,1e0,1e2,1e4])
         set(gca,'ColorScale','log')
         %clim(ax2,[0,CTopHIV])
         xlabel('Time (days)','FontSize',40,'FontWeight','Bold')
         ylabel('Depth (cm)','FontSize',40,'FontWeight','Bold')
         XT = get(gca,'XTick');
         set(gca,'FontSize',40)
         YT = get(gca,'YTick');
         set(gca,'FontSize',40)
         %ah = axes;
         hold on
         plot(ax,lineLE_X,lineLE_Y,'LineWidth',2,'Color',lineColor,'LineStyle','--')
         plot(ax,lineES_X,lineES_Y,'LineWidth',2,'Color',lineColor,'LineStyle','--')
         
      % else
      % 
      %    figure()
      %    title('Viral Load Distribution','FontSize',22,'FontWeight','Bold')
      %    surf(XX,YY,ZZ_2,'EdgeColor','none');
      %    view([0,90])
      %    set(gca,'Ydir','reverse')
      %    xlim([0,endT])
      %    ylim([0,h_E+h_S])
      %    colormap('pink')
      %    colorbar
      %    caxis([0,CTopHIV])
      %    xlabel('Time (hrs)','FontSize',18,'FontWeight','Bold')
      %    ylabel('Depth into Tissue (cm)','FontSize',18,'FontWeight','Bold')
      %    XT = get(gca,'XTick');
      %    set(gca,'FontSize',16)
      %    YT = get(gca,'YTick');
      %    set(gca,'FontSize',16)
      %    %ah = axes;
      %    hold on;
      % 
      % end
     
%      if ~isempty(timeINF)
%         %timeINFvec =  timeINF*ones(length(t),length(x));
%         plot3([timeINF timeINF],[0,0.3],[10^5,10^5],'--','LineWidth',2,'Color',[0.9,0.9,0.9]) % vertical line
%         %line(t,x,timeINF*ones(1,length(x)))
%         %surf(x,t,timeINFvec)
%         %colormap('gray')
%      end

     
end