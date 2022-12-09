function createLinePlots(t, Mesh, T_params, vir, chem, tar_cell, inf_cell, vir_b, tar_b, inf_b, settings)

    %% LOCAL - PLOTTING VIRAL LOAD VS. TIME
    if settings("showVvsTime") == 1
        
        f = figure();
        ax = axes('Parent',f);
        
        VL_avg = (trapz(vir(1:Mesh.numL,:),1) / Mesh.numL)';
        VT_avg = (trapz(vir(Mesh.numL+1:Mesh.numX,:),1) / Mesh.numT)';
          
        %plot(t/24,log10(VS_avg),'LineWidth',2)
        semilogy(ax, t/24, VL_avg, 'color', 'b', 'LineWidth',5)
        hold on
        semilogy(ax, t/24, VT_avg, 'color', 'r', 'LineWidth',5)
        %xlim([0 60])
        xlim([0, t(end)/24])
        ylim(ax, [10^-1 inf])
        %title('Viral Load in Stroma vs. Time','FontSize',22,'FontWeight','Bold')
        xlabel(ax, 'Time (days)','FontSize',32,'FontWeight','Bold')
        ylabel(ax, 'Viral Load (virions/ml)','FontSize',32,'FontWeight','Bold')
        legend(ax, {'Lumen', 'Tissue'})
        set(ax,'FontSize',28)
    end

    %% LOCAL - PLOTTING VIRAL LOAD VS R IN TISSUE @ DIFFERENT TIMES

    if settings("showVvsX") == 1
        %tVals = [1*60, 2*3600, 12*3600, 24*3600, 48*3600, 168*3600];
        tVals = [0*60, 6*60, 30*60, 2*60*60];
        tInds = zeros(1,length(tVals));
        legendLabels = cell(1,length(tVals));
    
        rLx = Mesh.x(Mesh.numL)*ones(5,1);
        rLy = linspace(0.001, max(max(vir)),5);

        rEx = Mesh.x(Mesh.numL+Mesh.numE)*ones(5,1);
        rEy = linspace(0.001, max(max(vir)),5);
    
        cm = colormap('turbo');
        cmInds = floor(linspace(1,length(cm),length(tVals)));
    
        for i = 1:length(tVals)
            tInds(i) = find(t*3600 >= tVals(i),1);
            legendLabels{i} = sprintf('%.1f hrs', tVals(i)/(3600));
        end
    
        f = figure();
        ax = axes('Parent',f);
        for i = 1:length(tInds)
            plot(ax, Mesh.x, vir(:,tInds(i)),'LineWidth',5)
            hold on
        end
        hold on
        plot(ax, rLx, rLy, 'Color',[0.7,0.7,0.7],'LineStyle','--','LineWidth',2)
        plot(ax, rEx, rEy, 'Color',[0.7,0.7,0.7],'LineStyle','--','LineWidth',2)
        set(ax,'FontSize',28)
        legend(ax,legendLabels)
        xlabel(ax,'Depth (cm)','FontWeight','Bold','FontSize',32)
        ylabel(ax,'Viral Load in Tissue (virions/cm^3)','FontWeight','Bold','FontSize',32)
        %ylim(ax,[1e-4,1e4])
        ylim(ax,[0, 1e3]);
        %xlim(ax,[0, Mesh.x(end)])
        xlim(ax,[0, 0.2])
    end
    
    %% LOCAL - PLOTTING CHEMOKINE VS. TIME

    if settings("showUvsTime") == 1
        f = figure();
        ax = axes('Parent',f);
        UT_avg = (trapz(chem,1) / Mesh.numT)';
        semilogy(ax, t/24,UT_avg, 'color', 'b', 'LineWidth',5)
        xlabel(ax, 'Time (days)','FontSize',32,'FontWeight','Bold')
        ylabel(ax, 'Chemokine Concentration (ml^{-1})','FontSize',32,'FontWeight','Bold')
        set(ax,'FontSize',28)
        ylim([-inf inf])
    %     xlim([0 60])
    %     ylim([10^-6 inf])
    end
    
    %% LOCAL - PLOTTING CELLS VS. TIME

    if settings("showCellsvsTime") == 1
        TT_avg = (trapz(tar_cell,1) / Mesh.numT)';
        IT_avg = (trapz(inf_cell,1) / Mesh.numT)';
   
        %format long
        %TT_avg(find(TT_avg<10^4,1))
        
        f = figure();
        ax = axes('Parent',f);
        plot(ax, t/24,[TT_avg,IT_avg],'LineWidth',5)
        xlim([0, t(end)/24])
        ylim(ax, [0, T_params.T0_T])
    %     ylim([0,max(max([TT_avg;IT_avg]))])
    %     xlim([0,60])
        %set(gca,'YScale','log');
        legend(ax, {'Target', 'Infected'})
        xlabel(ax, 'Time (days)','FontSize',32,'FontWeight','Bold')
        ylabel(ax, 'Cells Densities in Tissue (cells mg^{-1})','FontSize',32,'FontWeight','Bold')
        %title('Cell density vs. time','FontSize',22,'FontWeight','Bold')
        set(ax,'FontSize',28)
    end
    
    %% SYSTEMIC - PLOTTING V VS TIME

    if settings("showVvsTime_b") == 1
        f = figure();
        ax = axes('Parent',f);
        semilogy(ax, t/24, vir_b,'LineWidth',5)
        xlabel(ax, 'Time (days)','FontSize',32,'FontWeight','Bold')
        ylabel(ax, 'Viral Load in Blood (virions/ml)','FontSize',32,'FontWeight','Bold')
        set(ax,'FontSize',28)
        xlim([0, t(end)/24])
        max_vir_b = round(max(vir_b),-1);
        if max_vir_b == 0
            ylim([1e0, 1e1])
        else
            ylim([1e0, max_vir_b])
        end
    end

    %% SYSTEMIC - PLOTTING CELLS VS TIME

    if settings("showCellsvsTime_b") == 1
        f = figure();
        ax = axes('Parent',f);
        plot(ax, t/24,[tar_b',inf_b'],'LineWidth',5)
        xlabel(ax, 'Time (days)','FontSize',32,'FontWeight','Bold')
        ylabel(ax, 'Cell Densities in Blood (cells/ml)','FontSize',32,'FontWeight','Bold')
        set(ax,'FontSize',24)
        xlim([0, t(end)/24])
        legend(ax, {'Target', 'Infected'})

    end