% Modified and expanded from code written by Kent Ueno
% Linta Joseph; February 2025 
% 
% Description: 
%

addpath(genpath('QETLAB-0.9'));
maxNumCompThreads(16);

% Name of the file to store results 
filename = ['2025-03-04_EES_XZdisorder_6-6_spins_' ...
    'spin_op_midthird_u2=0.25-v2_J=32.76_h=6_g=0_v=0_integrable'];

% Xuniformfield=0.001h
J = 32.76; % natural dipolar coupling strength; krad/s 
h = 6.13; % natural disorder strength; krad/s

% g_vals = [0,0.02,0.004,0.02,0.06,0.14];%.*(1/h);
% v_vals = [0,0,0.16,0.01,0.01,0.01];%.*(1/J);

gr = 0;%(0:0.01:1.22).*(1/h);%0.02;%(0:0.02:0.5);%linspace(0,5,50).*(1/h);%(0:0.1:5).*(1/h);%(0:0.01:1.22).*(1/h);%0;(0:0.5:9).*(1/h);%0.1.*(1/h);%(0.01:1:50).*(1/h);%(0.1:0.2:4).*(1/h);%0.01:0.05:4;%0.1:.1:9;
vr = 0;%(0:0.05:0.5).*(1/J);%(0:0.01:0.2);%linspace(0,5,50).*(1/J);%(0:0.1:5).*(1/J);%(0:0.01:1).*(1/J);%0;(0:0.05:0.5).*(1/J);%0.1.*(1/J);;%(0:0.1:1);%(0:0.01:0.1).*(1/J);%0:.01:1;

vrl = length(vr);

cell2D = cell(length(gr),length(vr));
eigvals = cell(length(gr),length(vr));
EE = zeros(length(gr),length(vr));
EEbyEE_page = EE;

% number of spins
M = 6;
N = 6;

spectrum_third = round((2^(M+N))/3);

poolobj=parpool('local',16); % start at parallel pool of 16 cores on the local node
parfor m = 1:length(gr)
	for n = 1:vrl
         %u = 0.25/J;
         %u=0.24;
         u = sqrt(0.25-((J*vr(n))^2))*(1/J);

        [EE(m,n),eigvals{m,n},cell2D{m,n}] = Schmidt_ES_NMR(M,N,u,vr(n),gr(m),0,gr(m),J,h);
	    %cell2D{m,n} = Schmidt_ES_NMR(M,N,u,vr(n),gr(m),0,gr(m),J,h);

    end
end

delete(poolobj) % close the pool
clear poolobj vrl

save([filename,'_EES.mat'],'-v7.3')

m=length(gr);
n=length(vr);

% matrices to store fit values
fitmat = sparse(m,n);
distri_avg = zeros(m,n);
hist_vals_names = {'rowwise','x_eg','h_eg'};  
std_fns_norm_names = {'rowwise','P','SP','GOE','GUE'};

poolobj = parpool('local',16);
parfor i = 1:m
	for j = 1:n

	r_tot = [];						% total array of ratios
    ent_spectrum_midthird= [];
	% pick out the array of zets for an eigenvector, calculate ratios, then add to total ratio array r_tot
	for c = 1:length(cell2D{i,j})

	    zet = sort(cell2D{i,j}{1,c});						% sorted array of zetas
        
        if (length(vr)==1 && length(gr)==1) %for larger number of spins, save the ent spectrum DoS, not for the full phase space run
            ent_spectrum_midthird = horzcat(ent_spectrum_midthird,zet');
            [ent_dos,ent] = histcounts(ent_spectrum_midthird,300);%,'Normalization','probability');
	        ent = ent + (ent(2)-ent(1))/2;
	        ent(:,end) = [];
            ent_dos_hist{i,j} = [ent
                                 ent_dos];
        end

        if length(zet)<3
        else
	    r = sparse(1,(length(zet)-2));       				% ratios of consecutive spacings

	    % min/max ratio
	    for k = 2:(length(zet)-1)
	    	 r(k-1) = min([(zet(k+1)-zet(k)),(zet(k)-zet(k-1))])/max([(zet(k+1)-zet(k)),(zet(k)-zet(k-1))]);
        end

	    % store ratios from each eigenvector into total array
		r_tot = horzcat(r_tot,r);
        end
    end
    
    %Calculate average of the distribution -- 0.53 for GOE distribution
    distri_avg(i,j) = mean(r_tot);

	% determine histogram bin values and bin centers
	[rhor,r2] = histcounts(r_tot,300,'Normalization','pdf');
	r2 = r2 + (r2(2)-r2(1))/2;
	r2(:,end) = [];

	%rhor = rhor/2;
    
    x_eg = r2;
    h_eg = rhor;

    if (length(vr)==1 && length(gr)==1)
        hist_vals{i,j}= [x_eg
                         h_eg]; % saving the histogram; added 2025-02-06
    end

    % % fit functions (2 P(r))
	% ffit = 2*(1+r2).^(-2);    % Poisson
    % % ffit=(1/sum(ffit))*ffit;
	% ffit2 = 2*6.*r2./((1+r2).^(4));	% Semi-Poisson
    % % ffit2=(1/sum(ffit2))*ffit2;    
	% gfit = 2*((8/27)^(-1)).*(r2+r2.^2)./((1+r2+r2.^2).^2.5);  % GOE   
    % % gfit=(1/sum(gfit))*gfit;    
	% gfit2 = 2*((4*pi/81/sqrt(3))^(-1)).*((r2+r2.^2).^2)./((1+r2+r2.^2).^4); % GUE
    % % gfit2=(1/sum(gfit2))*gfit2;
    % 
    % if (length(vr)==1 && length(gr)==1)
    %     std_fns_norm{i,j} = [ffit
    %                         ffit2
    %                         gfit
    %                         gfit2]; %saving normalized standard funcitons 
    %                            % added 2025-02-07
    % end 
    % 
	% % fit to linear combination of SP and GOE
    % myfittype_lin = fittype('a*2*(6*r2/((1+r2)^4))+c*(((8/27)^(-1))*(r2+r2^2)/((1+r2+r2^2)^2.5))','dependent',{'rho_r'},'independent',{'r2'},'coefficients',{'a','c'});
    % [myfit_lin,gof] = fit(r2',rhor',myfittype_lin,'StartPoint',[0,0],'Lower',[0,0],'Upper',[1.5,1.5]);
    % 
    % coeff = coeffvalues(myfit_lin);
    % 
    % % fit to Poisson
    % myfittype_lin = fittype('a*2*(1+r2).^(-2) ','dependent',{'rho_r'},'independent',{'r2'},'coefficients',{'a'});
    % [myfit_lin_P,gof_P] = fit(r2',rhor',myfittype_lin,'StartPoint',[0],'Lower',[0],'Upper',[1.5]);
    % 
    % % fit to GUE
    % myfittype_lin = fittype('(a*2*(4*pi/81/sqrt(3))^(-1)).*((r2+r2.^2).^2)./((1+r2+r2.^2).^4) ','dependent',{'rho_r'},'independent',{'r2'},'coefficients',{'a'});
    % [myfit_lin_GUE,gof_GUE] = fit(r2',rhor',myfittype_lin,'StartPoint',[0],'Lower',[0],'Upper',[1.5]);
    % 
    % % fit to GOE
    % myfittype_lin = fittype('a*2*(((8/27)^(-1))*(r2+r2^2)/((1+r2+r2^2)^2.5)) ','dependent',{'rho_r'},'independent',{'r2'},'coefficients',{'a'});
    % [myfit_lin_GOE,gof_GOE] = fit(r2',rhor',myfittype_lin,'StartPoint',[0],'Lower',[0],'Upper',[1.5]);
    % 
    % % fit to Semi-Poisson
    % myfittype_lin = fittype('a*2*(6*r2/((1+r2)^4))','dependent',{'rho_r'},'independent',{'r2'},'coefficients',{'a'});
    % [myfit_lin_SP,gof_SP] = fit(r2',rhor',myfittype_lin,'StartPoint',[0],'Lower',[0],'Upper',[1.5]);
    % 
    % % calculate normalized fit to GOE
    % fitmat(i,j) = coeff(2)/(coeff(1)+coeff(2));%coeff(2)/(coeff(1)+coeff2));
    % fitmat_SP(i,j) = coeff(1)/(coeff(1)+coeff(2));%coeff(2)/(coeff(1)+coeff(2));
    % 
    % r_sq(i,j) = gof.rsquare;
    % r_sq_P(i,j) = gof_P.rsquare;
    % r_sq_SP(i,j) = gof_SP.rsquare;
    % r_sq_GUE(i,j) = gof_GUE.rsquare;
    % r_sq_GOE(i,j) = gof_GOE.rsquare;
    % 
    % JS_P(i,j) = DJS_P;
    % JS_SP(i,j) = DJS_SP;
    % JS_GOE(i,j) = DJS_GOE;
    % JS_GUE(i,j) = DJS_GUE;

    end
end

delete(poolobj);
clear poolobj;

save([filename,'_EES.mat'],'-v7.3') %moved to the end

if (length(vr)>1 && length(gr)>1) %only for full phase space runs
    % figure
    % subplot(2,3,1)
    % pcolor(vr,gr,fitmat)
    % xlabel('v')
    % ylabel('g')
    % cb = colorbar;
    % ylabel(cb,'coeff-GOE/(coeff-SP + coeff-GOE)')
    % %clim([0,1])
    % subplot(2,3,2)
    % pcolor(J*vr,gr*h,r_sq)
    % cb = colorbar;
    % xlabel('Jv')
    % ylabel('hg')
    % ylabel(cb,'rsquare')
    % %clim([0,1])
    % subplot(2,3,3)
    % pcolor(vr,gr,r_sq_P)
    % cb = colorbar;
    % xlabel('v')
    % ylabel('g')
    % ylabel(cb,'rsquare-P')
    % %clim([0,1])
    % subplot(2,3,4)
    % pcolor(vr,gr,r_sq_SP)
    % cb = colorbar;
    % xlabel('v')
    % ylabel('g')
    % ylabel(cb,'rsquare-SP')
    % %clim([0,1])
    % subplot(2,3,5)
    % pcolor(vr,gr,r_sq_GUE)
    % cb = colorbar;
    % xlabel('v')
    % ylabel('g')
    % ylabel(cb,'rsquare-GUE')
    % %clim([0,1])
    % subplot(2,3,6)
    % pcolor(vr,gr,r_sq_GOE)
    % cb = colorbar;
    % xlabel('v')
    % ylabel('g')
    % ylabel(cb,'rsquare-GOE')
    % %clim([0,1])
    % saveas(gcf,[filename,'_FITS_EES.fig'])

    % % imagesc(vr,gr,fitmat)
    % % set(gca,'YDir','reverse')
    % figure
    % subplot(2,2,1)
    % %title('XZ, 7-5,u=0.24, J=1, h=1') %u=0.24
    % pcolor(vr,gr,abs(JS_P))
    % cb = colorbar;
    % % clim([0,0.07])
    % xlabel('v')
    % ylabel('g')
    % ylabel(cb,'JS_P')
    % subplot(2,2,2)
    % pcolor(vr,gr,abs(JS_SP))
    % cb = colorbar;
    % % clim([0,0.07])
    % xlabel('v')
    % ylabel('g')
    % ylabel(cb,'JS_{SP}')
    % subplot(2,2,3)
    % pcolor(vr,gr,abs(JS_GOE))
    % cb = colorbar;
    % % clim([0,0.07])
    % xlabel('v')
    % ylabel('g')
    % ylabel(cb,'JS_{GOE}')
    % subplot(2,2,4)
    % pcolor(vr,gr,abs(JS_GUE))
    % cb = colorbar;
    % % clim([0,0.07])
    % xlabel('v')
    % ylabel('g')
    % ylabel(cb,'JS_{GUE}')
    % saveas(gcf,'2025-02-08_XZ disorder_7-5_spins_spin_op_midthird_u=pt25-vsq_J=1_h=1_JS_timetest.fig')
  
    figure
    pcolor(vr,gr,distri_avg)
    xlabel('v')
    ylabel('g')
    cb = colorbar;
    ylabel(cb,'avg of the distribution')
    saveas(gcf,[filename,'_AVG_EES.fig'])
    
    figure
    EE_page = M*log(2)-(1/2); %page value for equal cut
    figure
    pcolor(vr,gr,EE./EE_page)
    xlabel('v')
    ylabel('g')
    cb = colorbar;
    ylabel(cb,'EE_{avg}./EE_{page}')
    saveas(gcf,[filename,'_EE_avg.fig'])
else
    EE_page = M*log(2)-(1/2);
    EEbyEE_page = EE/EE_page;
end

if (length(vr)==1 && length(gr)==1) %not for full phase space runs

    % matrices to store fit values
    fitmat_eig = zeros(m,n);
    distri_avg_eig = zeros(m,n);
    
    poolobj = parpool('local',16);
    parfor i = 1:m
	    for j = 1:n
    
	    r_tot_eig = [];						% total array of ratios
        energy_spectrum=[];
    
	    % for given set of eigenvalues, calculate ratios, then add to total ratio array r_tot
	     
	    for c = 1:length(eigvals{i,j})

	        zet_eig = sort(eigvals{i,j}{1,c});
        				    % sorted array of zetas 
            energy_spectrum = horzcat(energy_spectrum,zet_eig');
            [energy_dos,energy] = histcounts(energy_spectrum,100);%,'Normalization','probability');
	        energy = energy + (energy(2)-energy(1))/2;
	        energy(:,end) = [];
            energy_dos_hist{i,j} = [energy
                                         energy_dos];
            zet_eig = zet_eig(spectrum_third:2*spectrum_third); %only mid third considered for statistics 
        
            if length(zet_eig)<3
            else
	            r_eig = sparse(1,(length(zet_eig)-2));       				% ratios of consecutive spacings
        
	            % min/max ratio
	            for k = 2:(length(zet_eig)-1)
	    	        r_eig(k-1) = min([(zet_eig(k+1)-zet_eig(k)),(zet_eig(k)-zet_eig(k-1))])/max([(zet_eig(k+1)-zet_eig(k)),(zet_eig(k)-zet_eig(k-1))]);
	    	        % r_eig(k-1) = (zet_eig(k+1)-zet_eig(k))/(zet_eig(k)-zet_eig(k-1)); %simple ratios, changed 2025-02-09
                end
        
	        % store ratios from each eigenvector into total array
		    r_tot_eig = horzcat(r_tot_eig,r_eig);
            end
        end
        
        %Calculate average of the distribution -- 0.53 for GOE distribution
        distri_avg_eig(i,j) = mean(r_tot_eig);
    
	    % determine histogram bin values and bin centers
	    [rhor_eig,r2_eig] = histcounts(r_tot_eig,300,'Normalization','pdf');
	    r2_eig = r2_eig + (r2_eig(2)-r2_eig(1))/2;
	    r2_eig(:,end) = [];
    
	    %rhor_eig = rhor_eig/2; %moved the factor of 2 to the defn of std
        %fns
        
        x_eg = r2_eig;
        h_eg = rhor_eig;
        
        %if (isscalar(vr) && isscalar(gr)) %redundant
            hist_vals_eig{i,j}= [x_eg
                                h_eg]; % saving the histogram for single runs; added 2025-02-06
        %end
    
        % % normalized fit functions
	    % ffit = 2*(1+r2_eig).^(-2); % Poisson
        % % ffit=(1/sum(ffit))*ffit;
	    % ffit2 = 2*6.*r2_eig./((1+r2_eig).^(4)); % Semi-Poisson
        % % ffit2=(1/sum(ffit2))*ffit2;    
	    % gfit = 2*((8/27)^(-1)).*(r2_eig+r2_eig.^2)./((1+r2_eig+r2_eig.^2).^2.5); % GOE  
        % % gfit=(1/sum(gfit))*gfit;    
	    % gfit2 = 2*((4*pi/81/sqrt(3))^(-1)).*((r2_eig+r2_eig.^2).^2)./((1+r2_eig+r2_eig.^2).^4); 
        % % gfit2=(1/sum(gfit2))*gfit2; % GUE    
        % 
        % %if (isscalar(vr) && isscalar(gr)) %redundant
        %     std_fns_norm_eig{i,j} = [ffit
        %                             ffit2
        %                             gfit
        %                             gfit2]; %saving normalized standard funcitons 
        %                                     % added 2025-02-07
        % %end
        % 
	    % % fit to linear combination of SP and GOE
        % myfittype_lin = fittype('a*(6*r2_eig/((1+r2_eig)^4))+c*(((8/27)^(-1))*(r2_eig+r2_eig^2)/((1+r2_eig+r2_eig^2)^2.5))','dependent',{'rho_r'},'independent',{'r2_eig'},'coefficients',{'a','c'});
        % [myfit_lin,gof_eig] = fit(r2_eig',rhor_eig',myfittype_lin,'StartPoint',[0,0],'Lower',[0,0],'Upper',[1.5,1.5]);
        % 
        % coeff_eig = coeffvalues(myfit_lin);
        % 
        % % fit to Poisson
        % myfittype_lin = fittype('a*(1+r2_eig).^(-2) ','dependent',{'rho_r'},'independent',{'r2_eig'},'coefficients',{'a'});
        % [myfit_lin_P,gof_P_eig] = fit(r2_eig',rhor_eig',myfittype_lin,'StartPoint',[0],'Lower',[0],'Upper',[1.5]);
        % 
        % % fit to GUE
        % myfittype_lin = fittype('(a*(4*pi/81/sqrt(3))^(-1)).*((r2_eig+r2_eig.^2).^2)./((1+r2_eig+r2_eig.^2).^4) ','dependent',{'rho_r'},'independent',{'r2_eig'},'coefficients',{'a'});
        % [myfit_lin_GUE,gof_GUE_eig] = fit(r2_eig',rhor_eig',myfittype_lin,'StartPoint',[0],'Lower',[0],'Upper',[1.5]);
        % 
        % % fit to GOE
        % myfittype_lin = fittype('a*(((8/27)^(-1))*(r2_eig+r2_eig^2)/((1+r2_eig+r2_eig^2)^2.5)) ','dependent',{'rho_r'},'independent',{'r2_eig'},'coefficients',{'a'});
        % [myfit_lin_GOE,gof_GOE_eig] = fit(r2_eig',rhor_eig',myfittype_lin,'StartPoint',[0],'Lower',[0],'Upper',[1.5]);
        % 
        % % fit to Semi-Poisson
        % myfittype_lin = fittype('a*(6*r2_eig/((1+r2_eig)^4))','dependent',{'rho_r'},'independent',{'r2_eig'},'coefficients',{'a'});
        % [myfit_lin_SP,gof_SP_eig] = fit(r2_eig',rhor_eig',myfittype_lin,'StartPoint',[0],'Lower',[0],'Upper',[1.5]);
        % 
        % % calculate normalized fit to GOE
        % fitmat_eig(i,j) = coeff_eig(2)/(coeff_eig(1)+coeff_eig(2));%coeff(2)/(coeff(1)+coeff2));
        % fitmat_SP_eig(i,j) = coeff_eig(1)/(coeff_eig(1)+coeff_eig(2));%coeff(2)/(coeff(1)+coeff(2));
        % 
        % r_sq_eig(i,j) = gof_eig.rsquare;
        % r_sq_P_eig(i,j) = gof_P_eig.rsquare;
        % r_sq_SP_eig(i,j) = gof_SP_eig.rsquare;
        % r_sq_GUE_eig(i,j) = gof_GUE_eig.rsquare;
        % r_sq_GOE_eig(i,j) = gof_GOE_eig.rsquare;
        % 
        end
    end

    delete(poolobj);
    clear poolobj;
end

save([filename,'_EES.mat'],'-v7.3')

% commenting these out because with 12 spins and no averaging, 
% the energy statistics is not good at all 

% if (length(vr)>1 && length(gr)>1)
%     figure
%     subplot(2,3,1)
%     pcolor(vr,gr,fitmat_eig)
%     xlabel('v')
%     ylabel('g')
%     cb = colorbar;
%     ylabel(cb,'coeff-GOE/(coeff-SP + coeff-GOE)')
%     %clim([0,1])
%     subplot(2,3,2)
%     pcolor(J*vr,gr*h,r_sq_eig)
%     cb = colorbar;
%     xlabel('Jv')
%     ylabel('hg')
%     ylabel(cb,'rsquare')
%     %clim([0,1])
%     subplot(2,3,3)
%     pcolor(vr,gr,r_sq_P_eig)
%     cb = colorbar;
%     xlabel('v')
%     ylabel('g')
%     ylabel(cb,'rsquare-P')
%     %clim([0,1])
%     subplot(2,3,4)
%     pcolor(vr,gr,r_sq_SP_eig)
%     cb = colorbar;
%     xlabel('v')
%     ylabel('g')
%     ylabel(cb,'rsquare-SP')
%     %clim([0,1])
%     subplot(2,3,5)
%     pcolor(vr,gr,r_sq_GUE_eig)
%     cb = colorbar;
%     xlabel('v')
%     ylabel('g')
%     ylabel(cb,'rsquare-GUE')
%     %clim([0,1])
%     subplot(2,3,6)
%     pcolor(vr,gr,r_sq_GOE_eig)
%     cb = colorbar;
%     xlabel('v')
%     ylabel('g')
%     ylabel(cb,'rsquare-GOE')
%     %clim([0,1])
%     saveas(gcf,[filename,'_FITS_energyspectrum.fig'])
% 
%     figure
%     pcolor(vr,gr,distri_avg_eig)
%     xlabel('v')
%     ylabel('g')
%     cb = colorbar;
%     ylabel(cb,'avg of the distribution')
%     saveas(gcf,[filename,'_AVG_energyspectrum.fig'])
% end
% 
