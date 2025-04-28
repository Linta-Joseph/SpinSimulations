% Modified and expanded from code written by Kent Ueno
% Linta Joseph; February 2025 
% 
% Description: Calculate the entanglement entropy for an initial product
% state for evolution under a linear spin chain Hamiltonian with controllable 
% interactions and disorder 

% Add the Qetlab package folder to Matlab path 
addpath(genpath('QETLAB-0.9'));
maxNumCompThreads(16);

% Number of spins, two subsystems of a linear chain
lp = 5;
rp = 5;

J = 32.76; % natural dipolar coupling strength; krad/s 
h = 6.13; % natural disorder strength; krad/s

%Experimentally controllable disorder and interaction strengths
gs = [0.02,0.16];%[0,0.02,0.02,0.06,0.14];
vs = [0.01,0.004];%[0,0,0.01,0.01,0.01];

% Name of the file to store results 
filename = ['2025-02-22_EE_6-6_50reals_XZdisorder_u2=0.25-v2_' ...
    'J=32.76_h=6.1_100reals_(0.5)t=1000_v=0.01,0.004_g=0.02,0.16'];

% Define the times and timesteps for time evolution
J_eff = 0.5; %J(sqrt(u2+v2))=0.5
tstart = 0.1/J_eff; %Jv is the correct interaction strength here 
tfin = 1000.1/J_eff; 
tinc=0.1/J_eff;
t = tstart:tinc:tfin;
tnum = length(t);

sym = 'xz';

reals = 50; % number of disorder realizations

% cells to store entanglement entropy results
EE = cell(length(gs),length(vs));

for m = 1:length(gs)   
    EE{m} = ESS_and_EE(lp,rp,vs(m),gs(m),tstart,tinc,tfin,reals,sym,tnum,J,h);
end

clear m n
clear i j k ans 

save([filename,'.mat'], '-v7.3');

% Page entropy
EE_page = (lp*log(2)-0.5);

% Plot EE vs time
for m = 1:length(EE)
    plot(J_eff*t,EE{m}./EE_page,'LineWidth',2, ...
            'DisplayName',['hg = ',num2str(h*gs(m)),', Jv = ',num2str(J*vs(m))])
    hold on
end

legend
xlabel('J_{eff}t')
ylabel('S(\rho_{A}(t))/S_{Page}')
set(gca,'YScale','linear')
set(gca,'XScale','log')
saveas(gcf,[filename,'.fig'])