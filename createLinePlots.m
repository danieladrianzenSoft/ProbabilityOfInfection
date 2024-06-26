function createLinePlots(t, Mesh, T_params, vir, drug, drugdp, chem, tar_cell, inf_cell, vir_b, drug_b, drugdp_b, tar_b, inf_b, drugType, settings)

    lineWidth = 12;
    if t(end)/24 > 18
        tend = 18;
    else
        tend = t(end)/24;
    end

    IC50 = 90;
    calcQ = @(c) 1./(1+(IC50./(c))); 
    %% LOCAL - PLOTTING VIRAL LOAD VS. TIME
    if settings("showVvsTime") == 1
        f = figure();
        f.Position = [383,363,793,627];
        ax = axes('Parent',f);
        
        cm = getCustomColormap(3);
        VL_avg = (trapz(vir(1:Mesh.numL,:),1) / Mesh.numL)';
        VT_avg = (trapz(vir(Mesh.numL+1:Mesh.numX,:),1) / Mesh.numT)';
          
        %plot(t/24,log10(VS_avg),'LineWidth',2)
        %plot(ax, t/24, VL_avg, 'LineWidth',10)
        semilogy(ax, t/24, VL_avg, 'color', cm(1,:), 'LineWidth',lineWidth)
        hold on
        %plot(ax, t/24, VT_avg, 'LineWidth',10)
        semilogy(ax, t/24, VT_avg, 'color', cm(2,:), 'LineWidth',lineWidth)
        semilogy(ax, t/24, vir_b,'LineWidth',lineWidth,'color',cm(3,:))
        %xlim([0 60])
        xlim([0, tend])
        %ylim(ax, [10^-1 max(max([VL_avg;VT_avg]))])
        ylim(ax, [1e-2,1e7])
        yticks([1e-2,1e0,1e2,1e4,1e6])
        %title('Viral Load in Stroma vs. Time','FontSize',22,'FontWeight','Bold')
        xlabel(ax, 'Time (days)','FontSize',40,'FontWeight','Bold')
        ylabel(ax, 'Conc (vir/ml)','FontSize',40,'FontWeight','Bold')
        title(ax, 'Viral Load','FontSize',40,'FontWeight','Bold')
        legend(ax, {'Lumen', 'Tissue','Blood'},"Location","northwest")
        legend boxoff
        set(ax,'FontSize',40)
    end

    %% LOCAL - PLOTTING VIRAL LOAD VS R IN TISSUE @ DIFFERENT TIMES

    if settings("showVvsX") == 1
        %tVals = [1*60, 2*3600, 12*3600, 24*3600, 48*3600, 168*3600];
        tVals = [0*60, 6*60, 30*60, 2*60*60, 10*60*60];
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
    
        cm = getCustomColormap(length(tInds),"colormap","purple-teal");
        f = figure();
        ax = axes('Parent',f);
        for i = 1:length(tInds)
            plot(ax, Mesh.x, vir(:,tInds(i)),'LineWidth',5,'color',cm(i,:))
            hold on
        end
        hold on
        plot(ax, rLx, rLy, 'Color',[0.7,0.7,0.7],'LineStyle','--','LineWidth',2)
        plot(ax, rEx, rEy, 'Color',[0.7,0.7,0.7],'LineStyle','--','LineWidth',2)
        set(ax,'FontSize',40)
        legend(ax,legendLabels)
        legend boxoff
        xlabel(ax,'Depth (cm)','FontWeight','Bold','FontSize',40)
        ylabel(ax,'Conc (virions/cm^3)','FontWeight','Bold','FontSize',40)
        title(ax,'Viral Load in Tissue','FontWeight','Bold','FontSize',40)
        %ylim(ax,[1e-4,1e4])
        %ylim(ax,[0, 1e3]);
        ylim(ax,[0, max(max(vir(:,tInds)))]);
        %ylim(ax,[0,inf])
        %xlim(ax,[0, Mesh.x(end)])
        xlim(ax,[0, 0.2])
    end
    
    %% LOCAL - PLOTTING CHEMOKINE VS. TIME

    if settings("showUvsTime") == 1
        f = figure();
        f.Position = [383,363,793,627];
        ax = axes('Parent',f);
        cm = getCustomColormap(1);
        UT_avg = (trapz(chem,1) / Mesh.numT)';
        semilogy(ax, t/24,UT_avg, 'color', cm(1,:), 'LineWidth',lineWidth)
        xlabel(ax, 'Time (days)','FontSize',40,'FontWeight','Bold')
        ylabel(ax, 'Conc (ml^{-1})','FontSize',40,'FontWeight','Bold')
        title(ax, 'Chemokine in Tissue','FontSize',40,'FontWeight','Bold')
        set(ax,'FontSize',40)
        ylim([-inf inf])
    %     xlim([0 60])
    %     ylim([10^-6 inf])
    end

    %% LOCAL - PLOTTING DRUG VS. TIME

    if settings("showDvsTime") == 1
        if (drugType == 0 || drugType == 1)
            DR_avg = [];
            DL_avg = (trapz(drug(1:Mesh.numL,:),1) / Mesh.numL)';
            DT_avg = (trapz(drug(Mesh.numL+1:Mesh.numX,:),1) / Mesh.numT)';
            DE_avg = (trapz(drug(Mesh.numL+1:Mesh.numL+Mesh.numE,:),1) / Mesh.numE)';
            DS_avg = (trapz(drug(Mesh.numL+Mesh.numE+1:Mesh.numX,:),1) / Mesh.numS)';
        else
            DR_avg = (trapz(drug(1:Mesh.numR,:),1) / Mesh.numR)';
            DL_avg = (trapz(drug(Mesh.numR+1:Mesh.numR+Mesh.numL,:),1) / Mesh.numL)';
            DT_avg = (trapz(drug(Mesh.numR+Mesh.numL+1:Mesh.numR+Mesh.numX,:),1) / Mesh.numT)';
            DE_avg = (trapz(drug(Mesh.numR+Mesh.numL+1:Mesh.numR+Mesh.numL+Mesh.numE,:),1) / Mesh.numE)';
            DS_avg = (trapz(drug(Mesh.numR+Mesh.numL+Mesh.numE+1:Mesh.numR+Mesh.numX,:),1) / Mesh.numS)';
        end

        % f = figure();
        % ax = axes('Parent',f);
        % cm = getCustomColormap(3);
        % 
        % plot(ax,t,DL_avg*(1e3/287.213),'LineWidth',10,'color',cm(1,:))
        % hold on
        % plot(ax,t,DE_avg*(1e3/287.213),'LineWidth',10,'color',cm(2,:))
        % plot(ax,t,DS_avg*(1e3/287.213),'LineWidth',10,'color',cm(3,:))
        % 
        % ylim([0,max(max([DL_avg;DE_avg;DS_avg]*(1e3/287.213)))])
        % %set(gca,'YScale','log');
        % legend('Lumen','Epithelium', 'Stroma')
        % xlabel('Time (hrs)','FontSize',40,'FontWeight','Bold')
        % %ylabel('TFV (ng ml^{-1})','FontSize',40,'FontWeight','Bold')
        % ylabel('TFV (fmol/mg)','FontSize',40,'FontWeight','Bold')
        % set(gca,'FontSize',40)
        % xlim([0,24])
        
        f = figure();
        f.Position = [383,363,793,627];
        ax = axes('Parent',f);
        cm = getCustomColormap(3);
        plot(t,DS_avg*(1e3/287.213),'LineWidth',lineWidth,'color',cm(1,:))

        %hold on
        %plot(t,drug_b*(1e3/287.213),'LineWidth',lineWidth,'color',cm(2,:))

        %ylim([0,max(CS_avg)])
        %set(gca,'YScale','log');
        %legend('Stroma')
        %legend('Stroma','Blood')
        xlabel('Time (hrs)','FontSize',40,'FontWeight','Bold')
        %ylabel('TFV in Stroma (ng ml^{-1})','FontSize',40,'FontWeight','Bold')
        ylabel('Conc (fmol/mg)','FontSize',40,'FontWeight','Bold')
        %title('TFV','FontSize',40,'FontWeight','Bold')
        title('TFV in Stroma','FontSize',40,'FontWeight','Bold')
        set(gca,'FontSize',40)
        xlim([0,24])
        %yticks([1e-2,1e0,1e2,1e4,1e6])
        %yticks([0,0.5e5,1e5,1.5e5,2e5,2.5e5])
        yticks([0,4e5,8e5,12e5])
        ylim([0,12e5])
    end

    if settings("showDdpvsTime") == 1
        %TFV-DP
        DdpE_avg = (trapz(drugdp(1:Mesh.numE,:),1) / Mesh.numE)';
        DdpS_avg = (trapz(drugdp(Mesh.numE+1:Mesh.numE+Mesh.numS,:),1) / Mesh.numS)';       
        f = figure();
        f.Position = [383,363,793,627];
        ax = axes('Parent',f);
        cm = getCustomColormap(2);
        %plot(ax,t,DdpS_avg,'LineWidth',10,'color',cm(2,:));
        plot(ax,t,DdpS_avg*(1e3/447.17),'LineWidth',lineWidth,'color',cm(1,:));
        %hold on
        %plot(ax,t,drugdp_b*(1e3/447.17),'LineWidth',lineWidth,'color',cm(2,:));
        xlabel('Time (hrs)','FontSize',40,'FontWeight','Bold')
        %ylabel('TFV-DP in Stroma (ng ml^{-1})','FontSize',40,'FontWeight','Bold')
        ylabel('Conc (fmol/mg)','FontSize',40,'FontWeight','Bold')
        %title('TFV-DP','FontSize',40,'FontWeight','Bold')
        title('TFV-DP in Stroma','FontSize',40,'FontWeight','Bold')
        %legend('Stroma','Blood')
        %ylim([0,2000])
        %ylim([0,15000])
        %yticks([1e-10,1e-5,1e0,1e5])
        %yticks([0,500,1000,1500,2000])
        %yticks([0,5000,10000,15000])

        set(gca,'FontSize',40)
        xlim([0,24])
    end
    
    %% LOCAL - PLOTTING CELLS VS. TIME

    if settings("showCellsvsTime") == 1
        TT_avg = (trapz(tar_cell,1) / Mesh.numT)';
        IT_avg = (trapz(inf_cell,1) / Mesh.numT)';
   
        %format long
        %TT_avg(find(TT_avg<10^4,1))
        cm = getCustomColormap(2);
        f = figure();
        f.Position = [383,363,793,627];
        ax = axes('Parent',f);
        plot(ax, t/24,TT_avg/1e4,'LineWidth',lineWidth,'color',cm(1,:))
        hold on
        plot(ax, t/24,IT_avg/1e4,'LineWidth',lineWidth,'color',cm(2,:))
        xlim([0, tend])
        ylim(ax, [0, T_params.T0_T]/1e4)
    %     ylim([0,max(max([TT_avg;IT_avg]))])
    %     xlim([0,60])
        %set(gca,'YScale','log');
        legend(ax, {'Target', 'Infected'})
        legend boxoff
        xlabel(ax, 'Time (days)','FontSize',40,'FontWeight','Bold')
        title(ax, 'Cells in Tissue','FontSize',40,'FontWeight','Bold')
        ylabel(ax, 'Density (x10^{4} cells/ml)','FontSize',40,'FontWeight','Bold')
        %title('Cell density vs. time','FontSize',22,'FontWeight','Bold')
        set(ax,'FontSize',40)
    end
    
    %% SYSTEMIC - PLOTTING V VS TIME

    if settings("showVvsTime_b") == 1
        f = figure();
        %f.Position = [383,363,793,627];
        f.Position = [383,550,744,440];
        ax = axes('Parent',f);
        cm = getCustomColormap(3);
        semilogy(ax, t/24, vir_b,'LineWidth',lineWidth,'color',cm(3,:))
        %plot(ax, t/24, vir_b,'LineWidth',10)
        xlabel(ax, 'Time (days)','FontSize',40,'FontWeight','Bold')
        ylabel(ax, 'Conc (vir/ml)','FontSize',40,'FontWeight','Bold')
        %title(ax, 'Viral Load in Blood','FontSize',40,'FontWeight','Bold')
        set(ax,'FontSize',40)
        xlim([0, 80])
        %ylim([1e-1,max(max(vir_b))])
        %ylim([1e-4,1e7])
        ylim([1e2,1e8])
        yticks([1e3,1e5,1e7])
    end

    %% SYSTEMIC - PLOTTING DRUG VS TIME

    if settings("showDvsTime") == 1
        conc = 1e-4:1e0:1e6;
        q = calcQ(conc);
        f = figure();
        f.Position = [383,363,793,627];
        ax = axes('Parent',f);
        cm = getCustomColormap(1);
        semilogx(ax, conc, q, 'color', cm(1,:), 'LineWidth',lineWidth)
        xlabel('Conc (fmol/mg)','FontSize',40,'FontWeight','Bold')
        %ylabel('TFV in Blood (ng/ml)','FontSize',40,'FontWeight','Bold')
        title('TFV-DP Dose Response','FontSize',40,'FontWeight','Bold')
        ylabel('q','FontSize',40,'FontWeight','Bold')
        set(gca,'FontSize',40)
        xlim([1e-2,1e6])
        xticks([1e-2,1e0,1e2,1e4,1e6])
        ylim([0,1])
    end

    if settings("showDvsTime_b") == 1
        f = figure();
        f.Position = [383,363,793,627];
        ax = axes('Parent',f);
        cm = getCustomColormap(1);
        plot(ax, t, drug_b*(1e3/287.213), 'color', cm(1,:), 'LineWidth',lineWidth)
        xlabel('Time (hrs)','FontSize',40,'FontWeight','Bold')
        %ylabel('TFV in Blood (ng/ml)','FontSize',40,'FontWeight','Bold')
        title('TFV in Blood','FontSize',40,'FontWeight','Bold')
        ylabel('Conc (fmol/mg)','FontSize',40,'FontWeight','Bold')
        set(gca,'FontSize',40)
        xlim([0,24])
        ylim([0,40])
    end
    if settings("showDdpvsTime_b") == 1
        % Drug: TFV-DP
        f = figure();
        f.Position = [383,363,793,627];
        ax = axes('Parent',f);
        cm = getCustomColormap(1);
        plot(ax, t, drugdp_b*(1e3/447.17), 'color', cm(1,:), 'LineWidth',lineWidth)
        xlabel('Time (hrs)','FontSize',40,'FontWeight','Bold')
        title('TFV-DP in Blood','FontSize',40,'FontWeight','Bold')
        ylabel('Conc (fmol/mg)','FontSize',40,'FontWeight','Bold')
        set(gca,'FontSize',40)
        xlim([0,24])
    end
    if settings("showqvsTime_b") == 1
        % Drug: TFV-DP
        f = figure();
        f.Position = [383,363,793,627];
        ax = axes('Parent',f);
        cm = getCustomColormap(1);
        plot(ax, t, calcQ(drugdp_b), 'color', cm(1,:), 'LineWidth',lineWidth)
        xlabel('Time (hrs)','FontSize',40,'FontWeight','Bold')
        title('Drug Efficacy q in Blood','FontSize',40,'FontWeight','Bold')
        ylabel('Conc (fmol/mg)','FontSize',40,'FontWeight','Bold')
        set(gca,'FontSize',40)
        xlim([0,24])
    end

    %% SYSTEMIC - PLOTTING CELLS VS TIME

    if settings("showCellsvsTime_b") == 1
        f = figure();
        f.Position = [383,363,793,627];
        ax = axes('Parent',f);
        cm = getCustomColormap(2);
        plot(ax, t/24,tar_b'/1e4,'LineWidth',lineWidth,'color',cm(1,:))
        hold on
        plot(ax, t/24,inf_b'/1e4,'LineWidth',lineWidth,'color',cm(2,:))
        xlabel(ax, 'Time (days)','FontSize',40,'FontWeight','Bold')
        ylabel(ax, 'Density (x10^{4} cells/ml)','FontSize',40,'FontWeight','Bold')
        title('Cells in Blood','FontSize',40,'FontWeight','Bold')
        set(ax,'FontSize',40)
        xlim([0, tend])
        legend(ax, {'Target', 'Infected'}, "location", "southwest")
        legend boxoff

    end