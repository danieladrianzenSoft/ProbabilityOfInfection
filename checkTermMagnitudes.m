function checkTermMagnitudes(t, sigmapr_term, r_term, stax_term)

figure()
plot(t, sigmapr_term, 'LineWidth',2)
xlabel('Time (hrs)', 'FontSize', 14, 'FontWeight', 'Bold')
ylabel('\sigma'' ', 'FontSize', 14, 'FontWeight', 'Bold')

figure()
plot(t,r_term, 'LineWidth', 2)
xlabel('Time (hrs)', 'FontSize', 14, 'FontWeight', 'Bold')
ylabel('r', 'FontSize', 14, 'FontWeight', 'Bold')

figure()
plot(t,stax_term, 'LineWidth', 2)
xlabel('Time (hrs)', 'FontSize', 14, 'FontWeight', 'Bold')
ylabel('S_{tax}', 'FontSize', 14, 'FontWeight', 'Bold')

end