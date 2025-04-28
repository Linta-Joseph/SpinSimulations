%% Check:
% filename, N, couplingscount, test value count, num_samples
% list of tau and p1, errors, sequences
% load dipolar couplings and error samples from file for consistency

clear 
close all  

filename = ['2024-03-08_GaussianSampling(0,1)_8spins_16couplings_' ...
    '_P1_1us_Tau_4us_Coupmax_5KHz_10_Overrotvals_100samples_logspacing_AVGOVERROT.mat'];

%Set parameters
sequence_list = {'WHH', 'MREV8', 'MREV16', 'BR24', 'CORY48', 'YXX24', 'YXX48'};%, 'SED48'};%, 'BR24-pi'};
N = 8; % num of spins 
Couplings_num = 16; % how many different coupling matrices to average over
couplingsCount = Couplings_num;%for compatibility with older code
Cycle_num = 12; 
tnum = 1:1:Cycle_num; %cycle indices 
dim = 2^N;
num_sample = 100;%INCREASE TO??
TestValueCount = 10;

P1_vals = 1*10^-6; %linspace(0,10E-6,TestValueCount);%logspace(-8,-5,TestValueCount);%[0,2E-6];%[0,1,1.5E-6,2E-6,2.5E-6,3E-6];
Tau_vals = 4*10^-6; %linspace(2,10,TestValueCount)*10^-6;%logspace(-6,-5,TestValueCount);%[2,4,8,15]*10^-6;%10^-6;%Jtau
Delta_vals = 0; %[0,logspace(-1,3,TestValueCount)*5];;%linspace(-5000,5000,TestValueCount);%[0,logspace(-1,2,TestValueCount)*5];%[10,100,1000];[0,logspace(-1,2,TestValueCount)*5];%[0,100,500];%Delta/J is the relevant parameter
Phasetrans_vals = 0; %linspace(-0.1,0.1,TestValueCount);%[0,logspace(-3,-1,TestValueCount)*5];%0.05;%0.05;%0.2;%linspace(-0.2,0.2,TestValueCount); % not wrt J?
Overrot_vals = [0,logspace(-3,-1,TestValueCount)*5]; %linspace(-0.1,0.1,TestValueCount);%0.1/3;%0.05;%0.1;%linspace(-0.1,0.1,TestValueCount); % not wrt J?
Couplingmax_vals = 2*pi*[5]*10^3; %2*pi*[0,logspace(-1,4,TestValueCount)*5];%;%[0,logspace(-1,4,TestValueCount)*5];%*10^3;%2*pi*5*10^3;[0.5,1,2,5,10]*10^3;
load('Delta_rand_vals.mat','Delta_rand_vals_100');
Overrot_rand_vals = Delta_rand_vals_100; %SEPARATE SAMPLING FOR OVERROT?

Sequence_num = length(sequence_list);
P1_num = length(P1_vals);
Tau_num = length(Tau_vals);
Delta_num = length(Delta_vals);
% PhaseTrans_num = length(Phasetrans_vals);
% Overrot_num = length(Overrot_vals);
Couplingmax_num = length(Couplingmax_vals);

%Initialize data structures for storage
Unitary_fidelity = zeros(Sequence_num, P1_num,Tau_num, Couplingmax_num, ...
    num_sample, Couplings_num);
Unitary_fidelity_faircomp = zeros(Sequence_num, P1_num, Tau_num, Couplingmax_num, ...
    Delta_num,num_sample,Couplings_num);

%initialize common variables and spin operators
initVars 
initCollectiveObs %dimensionless spin operator

% Initialize and construct dipolar Hamiltonians
Hdips = cell(Couplings_num, 1); %couplings_num - NO OF COUPLINGS TO BE AVGD OVER.
%Hdips STORES A DIFFERENT DIPOLAR HAMILTONIAN IN EACH OF ITS CELLS 

%%  Load the saved N by N coupling matrices (for consistency)
load(['coupling_matrices_Gaussian(0,1)_',num2str(N),'_',num2str(couplingsCount),'.mat'],'coupling_matrices') 

% Construct each Hamiltonian
for j = 1:Couplings_num    
    dip = coupling_matrices{j};%/(2*pi); %multiplied by 2pi while constructing the Hamiltonians
    Hdips{j} = getHdip(N, dim, x, y, z, dip); %FN GENERATES DIP HAMS WITH RANDOM COUPLINGS SPECIFIED BY THE SYMMETRIC COUPLING MATRIX dip
    clear dip          
end    

% Simulate evolution for each sequence
for i = 1:length(sequence_list) 
    i
    %% Iterate over different parameter values                    
    for j = 1:length(P1_vals) 
        pulse = P1_vals(j); 
        f1 = 1/4/pulse; 
        for k = 1:length(Tau_vals)
            tau = Tau_vals(k);  
            
            for l = 1:length(Couplingmax_vals)
                coupling = Couplingmax_vals(l)/3;%prefactor multiplying Hdip

                for m = 1:length(Overrot_vals) %CHECK!
                    for n = 1:num_sample

                        Delta = 0; %Delta_rand_vals_100(n)*Delta_vals(m);%Delta_vals(m);
                        overrot = Overrot_rand_vals(n)*Overrot_vals(m); %Overrot_vals(m)
                        phaseTrans = 0; %Phasetrans_vals(m);     
                        rotationError = 1+overrot;
                        sequenceName = sequence_list{i}; %choose a sequence

                        %% Configure Test Sequence
                        sequence = getSequence(sequenceName,X,Y); % X AND Y are collective spin observables,
                        % DEFINED IN initCollectiveObjs 

                        %% Set the pulses, taus and cycle times
                        Pulses = sequence.Pulses;
                        Taus = tau * sequence.Taus;      
                        tCyc = sum(Taus); 
                        
                        %% Compute phase transient direction for each pulse                        
                        for editPulse = 1:length(sequence.Pulses)
                            sequence.Phtr_array{editPulse} = ComputeTransients(sequence.Pulses{editPulse},X,Y);
                        end      

                        %% Iterate over the couplings
                        %---------------------------------------------------------------------------
                                        
                        %% Construct H-system for each dipolar coupling matrix
                        for c=1:Couplings_num  % c IS couplingsCount INDEX                          
                            Hdip = Hdips{c};
                            Hsys = Hdip*coupling + Z*Delta;       
                                        
                            %% obtain experimental unitary
                            % THE ONLY DIFF BETWEEN THE TWO IS ADDING THE Hint EVOLUTION DURING
                            % THE PULSE TO THE EXP UNITARY - THIS IS WHAT FINITE PULSE LENGTH
                            % MEANS HERE...
                            testUnitaries = getTestUnitaries(sequence,Hsys,pulse,tau,f1,dim,overrot,phaseTrans);                           %pulse HERE IS THE PULSE LENGTH
                            expUnitary = testUnitaries{1};
                            deltaUnitary = testUnitaries{2};
            
                            %% set theoretical or ideal Unitary
                            idealHam = 0;       
                                        
                            %% Calculate unitary fidelity : Tr(U_th*U_exp)
                            idealUnitary = expm(-1i*idealHam); %speye(dim);
                            Unitary_fidelity(i,j,k,l,m,n,c) = abs(trace(idealUnitary' * expUnitary)/2^N); % (consider exp unitary to be U^(1/M); M=sequence length/tau to compare sequences of different lengths )
                            Unitary_fidelity_faircomp(i,j,k,l,m,n,c) = abs(trace(idealUnitary' * (expUnitary^(1/round(tCyc/tau))))/2^N);% (consider exp unitary to be U^(1/M); M=sequence length/tau to compare sequences of different lengths )     
                            
                          end      %end of iteration over coupling                                                                            
                        
                    end      %end of iteration over sample   
                end %end of iteration over Delta values
            end %end of iteration over coupling_max
         end %end of iteration over tau vals
    end %end of iteration over p1 vals
end %end of iteration over sequence

%% Average over couplings and samples
Unitary_fidelity_couplingavg = mean(Unitary_fidelity, 7); %average over coupling 
Unitary_fidelity_faircomp_couplingavg = mean(Unitary_fidelity_faircomp, 7); 
Unitary_fidelity_couplingandsampleavg = mean(Unitary_fidelity_couplingavg, 6); %average over coupling 
Unitary_fidelity_faircomp_couplingandsampleavg = mean(Unitary_fidelity_faircomp_couplingavg, 6); 

%% save values 
save(filename,...
 'N', 'Cycle_num', 'Couplings_num', 'sequence_list', 'num_sample', 'coupling', 'Hdips',...
 'Delta_vals', 'Tau_vals', 'P1_vals', 'Phasetrans_vals', 'Overrot_vals', 'Couplingmax_vals',...
 'C_avg_72tau_sampleavg', 'Unitary_fidelity_faircomp_couplingandsampleavg',...
 'Delta_rand_vals_100');
