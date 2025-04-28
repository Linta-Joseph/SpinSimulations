% Linta Joseph, Feburary 2025 
% 
% Description: Plot 
% (1) entanglement spectrum statistics, 
% (2) eigenspectrum statistics,
% (3) energy density of states, 
% (4) etanglement energy density of states
% for a single phase space point (v,g) from saved histograms.

clear 
filename = ['2025-03-03_EES_XZdisorder_6-6_spins_spin_op_midthird_u2=0.25-v2' ...
    '_J=32.76_h=6_v=0.004_g=0.16_MBL_EES'];
load([filename,'.mat'])

% v and g indices
m = 1; %g index
n = 1; %v index

% entanglement spectrum histogram
x = hist_vals{m,n}(1,:);
y = hist_vals{m,n}(2,:);

% eigen spectrum histogram
x1 = hist_vals_eig{m,n}(1,:);
y1 = hist_vals_eig{m,n}(2,:);

% energy density of states histogram
x2 = energy_dos_hist{m,n}(1,:);
y2 = energy_dos_hist{m,n}(2,:);

% etanglement energy density of states histogram
x3 = ent_dos_hist{m,n}(1,:);
y3 = ent_dos_hist{m,n}(2,:);

% Standard functions (2 P(r) for r_tilde; 
% r-- ratios of consecutive spacings, r_tilde - min/max ratios)
r2 = 0:0.01:1;
f = 2*(1+r2).^(-2);    % Poisson
f2 = 2*6.*r2./((1+r2).^(4));	% Semi-Poisson
g = 2*((8/27)^(-1)).*(r2+r2.^2)./((1+r2+r2.^2).^2.5);  % GOE   
g2 = 2*((4*pi/81/sqrt(3))^(-1)).*((r2+r2.^2).^2)./((1+r2+r2.^2).^4); % GUE

figure(1)
scatter(x,y,15,'k','filled')
hold on 
plot(r2,f,'LineWidth',2)
plot(r2,f2,'LineWidth',2);%,r2,g2)
plot(r2,g,'LineWidth',2)
legend('data','P','SP','GOE');%,'GUE')
title(['v = ',num2str(vr(n)),', g = ',num2str(gr(m))]);
xlim([0,1])
ylim([0,2])
xlabel('$\tilde{r}_{ent}$','Interpreter','Latex')
ylabel('$P(\tilde{r}_{ent})$','Interpreter','Latex')
saveas(gcf,[filename,'_EES.fig'])

figure(2)
scatter(x1,y1,15,'k','filled')
hold on 
plot(r2,f,'LineWidth',2)
plot(r2,f2,'LineWidth',2);%,r2,g2)
plot(r2,g,'LineWidth',2)
xlim([0,1])
ylim([0,2])
legend('data','P','SP','GOE');%,'GUE')
title(['v = ',num2str(vr(n)),', g = ',num2str(gr(m))]);
xlabel('$\tilde{r}_{eig}$','Interpreter','Latex')
ylabel('$P(\tilde{r}_{eig})$','Interpreter','Latex')
saveas(gcf,[filename,'_Energyspectrum.fig'])

figure(3)
scatter(x2,y2./sum(y2),15,'k','filled')
xlabel('$E$','Interpreter','Latex')
ylabel('$\rho_{E}$','Interpreter','Latex')
title(['v = ',num2str(vr(n)),', g = ',num2str(gr(m))]);
saveas(gcf,[filename,'_Energy_DoS_scaled.fig'])

figure(4)
scatter(x3,y3./sum(y3),15,'k','filled')
xlabel('$\epsilon$','Interpreter','Latex')
ylabel('$\rho_{\epsilon}$','Interpreter','Latex')
title(['v = ',num2str(vr(n)),', g = ',num2str(gr(m))]);
saveas(gcf,[filename,'_Ent_DoS_scaled.fig'])
