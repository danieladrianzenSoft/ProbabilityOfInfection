function createBloodPlots(t,B) 


% Drug: TFV
figure()
plot(t,B(:,1), 'color', 'k', 'LineWidth',2)
xlabel('Time (hrs)','FontSize',20,'FontWeight','Bold')
ylabel('TFV in Blood (ng/ml)','FontSize',18,'FontWeight','Bold')
set(gca,'FontSize',18)
xlim([0,24])
% Drug: TFV-DP
figure()
plot(t,B(:,2), 'color', 'k', 'LineWidth',2)
xlabel('Time (hrs)','FontSize',20,'FontWeight','Bold')
ylabel('TFV-DP in Blood (ng/ml)','FontSize',18,'FontWeight','Bold')
set(gca,'FontSize',18)
xlim([0,24])

% Cells:
figure()
plot(t/24,B(:,3), 'color', 'k', 'LineWidth',2)
xlabel('Time (days)','FontSize',20,'FontWeight','Bold')
%ylabel('Target Cells','FontSize',16,'FontWeight','Bold')
hold on
plot(t/24,B(:,4), 'color', 'r', 'LineWidth',2)
xlabel('Time (days)','FontSize',20,'FontWeight','Bold')
%ylabel('Latent Cells','FontSize',16,'FontWeight','Bold')
plot(t/24,B(:,5), 'color', 'b', 'LineWidth',2)
xlabel('Time (days)','FontSize',20,'FontWeight','Bold')
ylabel('Cell Densities in Blood (cells/ml)','FontSize',18,'FontWeight','Bold')
legend('Target Cells', 'Latent Cells', 'Infected Cells')
set(gca,'FontSize',18)


% Virus:
figure()
%plot(t/24,log10(B(:,6)), 'color', 'r', 'LineWidth',2)
semilogy(t/24,B(:,6), 'color', 'r', 'LineWidth',2)
xlabel('Time (days)','FontSize',20,'FontWeight','Bold')
ylabel('Viral Load in Blood (virions/ml)','FontSize',18,'FontWeight','Bold')
set(gca,'FontSize',18)
ylim([10^1, 10^5])
xlim([0,60])


% figure()
% plot(t,B(:,6), 'color', 'b', 'LineWidth',2)
% xlabel('Time (hrs)','FontSize',20,'FontWeight','Bold')
% ylabel('TFV-DP in Blood (ng/ml)','FontSize',18,'FontWeight','Bold')
% set(gca,'FontSize',18)
% xlim([0,24])



end