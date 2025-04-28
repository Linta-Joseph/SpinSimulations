clear 
filename = ['2025-03-03_EES_XZdisorder_6-6_spins_spin_op_midthird_u2=0.25-v2' ...
    '_J=32.76_h=6_v=0.004_g=0.16_MBL_EES'];
load([filename,'.mat'])

% v and g indices
m = 1; %g index
n = 1; %v index

x = hist_vals{m,n}(1,:);
y = hist_vals{m,n}(2,:);

x1 = hist_vals_eig{m,n}(1,:);
y1 = hist_vals_eig{m,n}(2,:);

x2 = energy_dos_hist{m,n}(1,:);
y2 = energy_dos_hist{m,n}(2,:);
%y2_prime = histogram()

x3 = ent_dos_hist{m,n}(1,:);
y3 = ent_dos_hist{m,n}(2,:);

% f = std_fns_norm{m,n}(1,:);
% f2 = std_fns_norm{m,n}(2,:);
% g = std_fns_norm{m,n}(3,:);
% g2 = std_fns_norm{m,n}(4,:);

r2 = 0:0.01:1;
% fit functions (2 P(r))
	f = 2*(1+r2).^(-2);    % Poisson
    % ffit=(1/sum(ffit))*ffit;
	f2 = 2*6.*r2./((1+r2).^(4));	% Semi-Poisson
    % ffit2=(1/sum(ffit2))*ffit2;    
	g = 2*((8/27)^(-1)).*(r2+r2.^2)./((1+r2+r2.^2).^2.5);  % GOE   
    % gfit=(1/sum(gfit))*gfit;    
	g2 = 2*((4*pi/81/sqrt(3))^(-1)).*((r2+r2.^2).^2)./((1+r2+r2.^2).^4); % GUE
    % gfit2=(1/sum(gfit2))*gfit2;


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

% annotation('textbox', [0.5, 0.2, 0.1, 0.1], 'String', '(b)',...
%     'FontSize',30,'BackgroundColor','None','EdgeColor','None',...
%     'Position',[-0.0074,0.968,0.057,0.044],'Units','normalized',...
%     'FontName','Times');
% 
% set(findall(gca,'-property','FontSize'),'FontSize',30,'FontName','Times');%,'FontWeight','Bold');
% legend('FontSize',25,'FontName','Times','EdgeColor','None',...
%     'Color','None','NumColumns',1)
% box on
% grid off
% 
% set(gca,'LineWidth',1)
% 
% f = gcf;
% f.Units = 'centimeters';
% %set(gcf,'Position',[0,0,18,12]) %3,43,838,593 %
% set(gcf,'Position',[0,0,24,18]) %3,43,838,593 %
% %set(gcf,'Position',[0,0,27,21]) %3,43,838,593 % for mqc figs with 2 subplots
% set(gcf,'paperposition',[0 0 4 3],'papersize', [8.5 11])
% %set(gcf,'OuterPosition',[-0.2117,-0.2117,24.42,20.45]) %3,43,838,593 %
% 
% %exportgraphics(gcf,'u=norm_EE.pdf')
