% Calculate the MQC width parameter and plot w.r.t time
%for different interaction and disorder strengths
clear 

tic

type = 'ZZ'; %type of otoc
sym = 'XZ'; %symmtery of the Hamiltonian
M = 4;
N = 4;
filename = ['2025-02-18_ZZ_OTOC_8spins_XZdisorder_J=33_h=6_Jt=1E5_vary_g'];


J = 32.76; % krad/s, natural interaction strength
hi = 6.13; %krad/s, natural disorder strength

% 4 g and v combinations
g_vals = linspace(0, 0.2, 10);
v_vals = 0.005;
dim = 2^(N+M); % Dimension of combined Hilbert space
coupling = 1; %NN coupling

phi_vals = linspace(-pi, pi, 128);
cohs = -64:1:63;

J_eff = J;

tstart1 = 1/J_eff; %Jv is the correct interaction strength here 
tfin1 = 10/J_eff; 
tinc1 = 1/J_eff;

tstart2 = 11/J_eff; %Jv is the correct interaction strength here 
tfin2 = 101/J_eff; 
tinc2 = 10/J_eff;

tstart3= 111/J_eff; %Jv is the correct interaction strength here 
tfin3 = 1011/J_eff; 
tinc3 = 100/J_eff;

tstart4= 1111/J_eff; %Jv is the correct interaction strength here 
tfin4 = 11111/J_eff; 
tinc4 = 1000/J_eff;

tstart5= 12111/J_eff; %Jv is the correct interaction strength here 
tfin5 = 112111/J_eff; 
tinc5 = 10000/J_eff;

t = horzcat(tstart1:tinc1:tfin1, tstart2:tinc2:tfin2,...
    tstart3:tinc3:tfin3, tstart4:tinc4:tfin4,...
    tstart5:tinc5:tfin5);
tnum = length(t);

n_reals = 100; %number of disorder realizations

n_g = length(g_vals);
n_phi = length(phi_vals);

% Pauli matrices
z=sparse([1 0; 0 -1]);
x=sparse( [ 0 1;1 0]); 
y=1i*sparse([0 -1;1 0]);
id=speye(2);

% Raising and lowering operators
p=sparse([0 1;0 0]);
m=sparse([0 0; 1 0]);
Sz=sparse(0.5.*[1 0; 0 -1]); Sx=sparse(0.5.*[ 0 1;1 0]); Sy=1i*sparse(0.5.*[0 -1;1 0]);

% Coupling matrix
% Nearest neighbor couplings
a=diag(ones(N+M-1,1),1) + diag(ones(N+M-1,1),-1);

% Z_collective is the collective nuclear spin polarization
Z_collective=sparse(dim, dim);
for k=1:N+M
    Z_collective = Z_collective + mykron(speye(2^(k-1)), Sz, speye(2^(N+M-k)));
end
X_collective=sparse(dim, dim);
for k=1:N+M
    X_collective = X_collective + mykron(speye(2^(k-1)), Sx, speye(2^(N+M-k)));
end

if (type == 'XZ')
    W_init = X_collective;
elseif (type == 'ZZ')
    W_init = Z_collective;
end

% Matrices to store widths
widths = zeros(tnum, n_reals);
widths_avg = zeros(tnum, n_g);

V = cell(1,length(phi_vals));
for phi_ind = 1:length(phi_vals)
    phi = phi_vals(phi_ind);
	V{phi_ind} = expm(1i * phi * Z_collective); % Phase shift operator
end

for i = 1:length(g_vals) %loop for every pair (g,v)
    i
    g = g_vals(i); %X disorder strength
    v = v_vals; %interaction strength
    u = sqrt(0.25 - (J*v)^2) .* 1/J;
    
    %Build up the Hamiltonian
    H=sparse(dim, dim);  
    
    for k=1:(N+M)
       for l=k+1:(N+M)
	    H=H+a(k,l)*J*(((u+v)/2)*mykron(speye(2^(k-1)),Sx,speye(2^(l-k-1)),Sx,speye(2^(N+M-l)))+...
	       ((v-u)/2)*mykron(speye(2^(k-1)),Sy,speye(2^(l-k-1)),Sy,speye(2^(N+M-l)))-...
	        v*mykron(speye(2^(k-1)),Sz,speye(2^(l-k-1)),Sz,speye(2^(N+M-l)))) ;

	     end
     end
     
     delete(gcp('nocreate'))
     poolobj=parpool('local',16); % start at parallel pool of 16 cores on the local node

     parfor reals_ind = 1:n_reals
        reals_ind
        Htot = sparse(dim, dim);

        X = sparse(dim, dim);
        Y = sparse(dim, dim);
        Z = sparse(dim, dim);   
        
        % build up disorder term
        for k=1:(N+M)           
            gxi = hi* g*(-1+2*rand);
            X = X + gxi*mykron(speye(2^(k-1)),Sx,speye(2^(N+M-k)));
            gyi = 0;%hi *g*(-1+2*rand);
            Y = Y + gyi*mykron(speye(2^(k-1)),Sy,speye(2^(N+M-k)));
            gzi = hi* g*(-1+2*rand);
            Z = Z + gzi*mykron(speye(2^(k-1)),Sz,speye(2^(N+M-k)));
        end

        % total Hamiltonian
        if sym == 'Z'
            Htot = H + Z;
        elseif sym == 'XZ'
            Htot = H + X + Z;
        elseif sym == 'XYZ'
            Htot = H + X + Y + Z;
        end
         
		% define Unitary operators for first and subsequent time steps
        % evolve from zero or to initial time if not from zero
        if tstart1 == 0

            % define unitary for time evolution
            U1=sparse(expm(-1i*Htot*tinc1));
            U2=sparse(expm(-1i*Htot*tinc2));
            U3=sparse(expm(-1i*Htot*tinc3));
            U4=sparse(expm(-1i*Htot*tinc4));
            U5=sparse(expm(-1i*Htot*tinc4));
            % initial state is a product state
            W = W_init;
      
         else
        
            % define unitaries for initial and subsequent time evolutions
            U0=sparse(expm(-1i*Htot*tstart1));
            
            % define unitary for subsequent time steps
            U1=sparse(expm(-1i*Htot*tinc1));  
            U2=sparse(expm(-1i*Htot*tinc2)); 
            U3=sparse(expm(-1i*Htot*tinc3)); 
            U4=sparse(expm(-1i*Htot*tinc4));
            U5=sparse(expm(-1i*Htot*tinc5));
        
            % evolve to initial time
            %check_unitarity = ctranspose(U)*U;
		    %W = U * Z_collective * U'; % Time evolve operator forward by time 
            W = U0 * W_init * U0'; %from original code
        
         end
    
         for t_index = 2:tnum    
            
            if(t(t_index)<=tfin1)
                U = U1;
            elseif(t(t_index)<=tfin2)
                U = U2;
            elseif(t(t_index)<=tfin3)
                U = U3;
            elseif(t(t_index)<=tfin4)
                U = U4;
            elseif(t(t_index)<=tfin5)
                U = U5;
            end

            C_phi = zeros(1, n_phi);

	        for phi_ind = 1:length(phi_vals)
	            C_phi(phi_ind) ... 
                    = trace(W' * (V{phi_ind})' * W * V{phi_ind}); % Calculate OTOC for given t and phi
            end

			C_phi_tot = [C_phi, C_phi(2:end), ...
                       C_phi(2:end)]; %WHAT'S THIS? WHY MAKE 3 COPIES..? 

			C_n_temp = fftshift(fft(C_phi_tot)); % Fourier transform OTOC data 
                     
			%% Peaks of each fourier transform, use abs bc the data has extremely small
            % % complex part % IS THIS STILL TRUE??
            [pks, locs] = findpeaks(abs(C_n_temp)); %pks and locs are =itermax X n_g X n_h X n_v cells
            locs1 = (locs - 192)/3;% 65th is 0 coherence
            pks = pks./sum(pks);
                    
             sum_temp = 0;
             for k = 1:length(pks)
                pk_num = locs1(k);
                sum_temp = sum_temp + pks(k) * pk_num^2; 
              end
              widths(t_index,reals_ind) = 2*(sum_temp);
     
              % evolve operator to the next timestep 
              W = U * W * U';
         end
     end
     clear U U0 Htot W
     delete(poolobj) % close the pool
     clear poolobj 

     clear H Htot %build a new Hamiltonian for each v
	 widths_avg(:,i) = mean(widths,2); % average over realizations
end

save([filename,'.mat'])

figure
for i=1:n_g
    plot(J*0.01*t,widths_avg(:,i),...
        'LineWidth',2,'DisplayName', ...
        ['g, h = ',num2str(g_vals(i)),', ' ...
        'v = ',num2str(v_vals(i))])
    hold on
end

xlabel('J_{eff}t')
ylabel('MQC widths (ZZ OTOC)')
set(gca,'YScale','linear')
set(gca,'XScale','log')
saveas(gcf,[filename,'.fig'])

toc