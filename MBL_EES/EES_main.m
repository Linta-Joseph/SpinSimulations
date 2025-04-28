% Modified and expanded from code written by Kent Ueno
% Linta Joseph; February 2025 
% 
% Description: Calculate the entanglement spectrum and eigenspectrum 
% statistics of a linear spin chain Hamiltonian with controllable 
% interactions and disorder 

% Add the Qetlab package folder to Matlab path 
addpath(genpath('QETLAB-0.9'));
maxNumCompThreads(16);

%% Set parameters 

% Name of the file to store results 
filename = ['2025-03-04_EES_XZdisorder_6-6_spins_' ...
    'spin_op_midthird_u2=0.25-v2_J=32.76_h=6_g=0_v=0_integrable'];

% Xuniformfield=0.001h
J = 32.76; % natural dipolar coupling strength; krad/s 
h = 6.13; % natural disorder strength; krad/s

%Experimentally controllable disorder and interaction strengths
% Run a full grid of values -- arrays for each gr and vr - for 12 spins OR
% Run a single value for each - for 14 spins
gr = (0:0.01:1.22).*(1/h);
vr = (0:0.05:0.5).*(1/J);

% number of spins - 2 subsystems of a linear chain
M = 6;
N = 6;

%%

vrl = length(vr);

%initialize cell arrays to store results 
cell2D = cell(length(gr),length(vr));
eigvals = cell(length(gr),length(vr));
EE = zeros(length(gr),length(vr));
EEbyEE_page = EE;

poolobj=parpool('local',16); % start at parallel pool of 16 cores on the local node
parfor m = 1:length(gr)
	for n = 1:vrl
         % set u with or without a chosen normlaization condition
         u = 0.25/J;
         %u=0.24;
         %u = sqrt(0.25-((J*vr(n))^2))*(1/J);

        [EE(m,n),eigvals{m,n},cell2D{m,n}] = Schmidt_ES_NMR(M,N,u,vr(n),gr(m),0,gr(m),J,h);

    end
end

delete(poolobj) % close the pool
clear poolobj vrl

m=length(gr);
n=length(vr);

% matrix to store results
distri_avg = zeros(m,n);

hist_vals_names = {'rowwise','x_eg','h_eg'};  
std_fns_norm_names = {'rowwise','P','SP','GOE','GUE'};

%% Calculations for the ENTANGLEMENT SPECTRUM 
poolobj = parpool('local',16);
parfor i = 1:m
	for j = 1:n

	    r_tot = [];						% full array of spacing ratios
        ent_spectrum_midthird= [];
	    % pick out the array of zets for an eigenvector, calculate ratios, 
        % then add to total ratio array r_tot
	    for c = 1:length(cell2D{i,j})
    
	        zet = sort(cell2D{i,j}{1,c});						% sorted array of zetas
                
            %for larger number of spins, only one point in phase space (v,g)
            % is run at a time, If so, save the ent spectrum DoS, not for the full phase space run
            if (length(vr)==1 && length(gr)==1) 
                ent_spectrum_midthird = horzcat(ent_spectrum_midthird,zet');
                [ent_dos,ent] = histcounts(ent_spectrum_midthird,300);
	            ent = ent + (ent(2)-ent(1))/2;
	            ent(:,end) = [];
                ent_dos_hist{i,j} = [ent
                                     ent_dos];
            end
    
            if length(zet)<3
            else
            % initialize array to store spacing ratios    
	        r = sparse(1,(length(zet)-2)); 
    
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
    
	    % determine entanglement spectrum histogram bin values and bin centers
	    [rhor,r2] = histcounts(r_tot,300,'Normalization','pdf');
	    r2 = r2 + (r2(2)-r2(1))/2;
	    r2(:,end) = [];
    
	    %rhor = rhor/2;
        
        x_eg = r2;
        h_eg = rhor;
        
        % saving the histogram; added 2025-02-06
        if (length(vr)==1 && length(gr)==1)
            %entanglement spectrum spacing ratio statistics histogram
            hist_vals{i,j}= [x_eg
                             h_eg]; 
        end

    end
end

delete(poolobj);
clear poolobj;

%% Entanglement Spectrum Plots
if (length(vr)>1 && length(gr)>1) %only for full phase space runs
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

%% Calculations for the EIGENSPECTRUM 
% Done only for single phasespace point runs, for M=N=7

if (length(vr)==1 && length(gr)==1) %not for full phase space runs

    % matrices to store fit values
    fitmat_eig = zeros(m,n);
    distri_avg_eig = zeros(m,n);
    
    poolobj = parpool('local',16);
    parfor i = 1:m
	    for j = 1:n
    
	    r_tot_eig = [];	% total array of ratios
        energy_spectrum=[];
    
	    % for given set of eigenvalues, calculate ratios, 
        % then add to total ratio array r_tot
	     
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
        
        x_eg = r2_eig;
        h_eg = rhor_eig;
        
        %eigenvalue spacing ratio statistics histogram
        hist_vals_eig{i,j}= [x_eg
                                h_eg]; 
    
        end
    end

    delete(poolobj);
    clear poolobj;
end

save([filename,'_EES.mat'],'-v7.3')
