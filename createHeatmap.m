function createHeatmap(t,x,DrugTissue,HIVTissue,h_E,h_S)

    endT = 120; %hours
    [XX,YY]=meshgrid(t,x);
    ZZ_2 = real(log10(HIVTissue));
    CTopHIV = floor(max(max(ZZ_2))) + floor( max(max(ZZ_2)-floor(max(max(ZZ_2))))/0.5) * 0.5;

      if max(max(DrugTissue)) > 0

         figure()
         tiledlayout(2,1)

         ZZ_1 = log10(DrugTissue);     
         CTopDrug = floor(max(max(ZZ_1))) + floor( max(max(ZZ_1)-floor(max(max(ZZ_1))))/0.5) * 0.5;
         title('TFV-DP Distribution','FontSize',22,'FontWeight','Bold')
         ax1 = nexttile;
         %contourf(XX,YY,ZZ_1,'LineColor','none')
         surf(XX,YY,ZZ_1,'EdgeColor','none')
         view([0,90])
         set(gca,'Ydir','reverse')
         xlim([0,endT])
         ylim([0,h_E+h_S])
         colormap(ax1,'bone')
         colorbar
         caxis([0,CTopDrug])
         xlabel('Time (hrs)','FontSize',18,'FontWeight','Bold')
         ylabel('Depth into Tissue (cm)','FontSize',18,'FontWeight','Bold')
         XT = get(gca,'XTick');
         set(gca,'FontSize',16)
         YT = get(gca,'YTick');
         set(gca,'FontSize',16)

         ax2 = nexttile;
         %contourf(XX,YY,ZZ_2,'LineColor','none')
         title('Viral Load Distribution','FontSize',22,'FontWeight','Bold')
         surf(XX,YY,ZZ_2,'EdgeColor','none')
         view([0,90])
         set(gca,'Ydir','reverse')
         xlim([0,endT])
         ylim([0,h_E+h_S])
         colormap(ax2,'pink')
         colorbar
         caxis([0,CTopHIV])
         xlabel('Time (hrs)','FontSize',18,'FontWeight','Bold')
         ylabel('Depth into Tissue (cm)','FontSize',18,'FontWeight','Bold')
         XT = get(gca,'XTick');
         set(gca,'FontSize',16)
         YT = get(gca,'YTick');
         set(gca,'FontSize',16)
         %ah = axes;
         hold on;
         
      else
          
         figure()
         title('Viral Load Distribution','FontSize',22,'FontWeight','Bold')
         surf(XX,YY,ZZ_2,'EdgeColor','none');
         view([0,90])
         set(gca,'Ydir','reverse')
         xlim([0,endT])
         ylim([0,h_E+h_S])
         colormap('pink')
         colorbar
         caxis([0,CTopHIV])
         xlabel('Time (hrs)','FontSize',18,'FontWeight','Bold')
         ylabel('Depth into Tissue (cm)','FontSize',18,'FontWeight','Bold')
         XT = get(gca,'XTick');
         set(gca,'FontSize',16)
         YT = get(gca,'YTick');
         set(gca,'FontSize',16)
         %ah = axes;
         hold on;
          
      end
     
%      if ~isempty(timeINF)
%         %timeINFvec =  timeINF*ones(length(t),length(x));
%         plot3([timeINF timeINF],[0,0.3],[10^5,10^5],'--','LineWidth',2,'Color',[0.9,0.9,0.9]) % vertical line
%         %line(t,x,timeINF*ones(1,length(x)))
%         %surf(x,t,timeINFvec)
%         %colormap('gray')
%      end

     
end