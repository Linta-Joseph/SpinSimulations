%% sequenceAnalysis.m
% Wynter Alford
% January 2022
%
% For a given pulse sequence, does the following:
%    - Calculates the a specified number of terms of the magnus
%    expansion for that sequence using the recursive method
%
%    - Calculates the fidelity of that sequence against H=0, H=H0,
%    and H=H0+H1+... up to the highest order term chosen
%
%    - Calculates fidelity when using experimental or instantaneous pulses
%
%    - Can analyze sequence behaviour and robustness with respect to 
%    offsets, tau spacing, pulse length, dipolar coupling strength,
%    overrotations, and phase transients


% Added comments + modified  - Linta Joseph 
% March 2022

clear 

%% Control Parameters (change these)

global dim Pulses knownOmegas knownQs knownPs 

% Valid sequences are: WHH, MREV8, CORY48, YXX48, YXX24, YXX24S, AZ48,
%    AS10M, AS10E, WHHWHH, ML10, WPW9-2Cycle, MG8
sequence_list = {'WHH', 'MREV8', 'MREV16', 'BR24','CORY48', 'YXX24', 'YXX48'};%,...
                 %'SED36','SED36err','SED48'};

testVarName = 'Tau'; % parameter to test over 
%CHOOSE FROM (Delta, Tau, Coupling, coupling_lo, phaseTrans,% delta_lo, 'overrot_hi')

varMaxMod = 1; % factor by which to multiply the MAX TEST VAR

mode = 'max'; %can be 'max' to get up to maxTerm terms,clear
%or 'time' to compute for a certain amount of time %???
maxTerm = 4; % highest Magnus series term to compute 
computationTime = 30; % once elapsed time reaches this many seconds, no further terms are computed for this loop.

testValueCount = 20; % HOW MANY DIFF VALS OF THE TEST VARIABLE
couplingsCount = 4; % how many different coupling matrices to average over

N = 4; % NUM OF SPINS
pulse = 0;%1.4e-6; % PULSE LENGTH (CLOSE TO EXPERIMENTAL)
tau = 7.4e-6;  % DELAY B/W PULSES IN S
coupling = 2*pi*420; %Hz
Delta = 2*pi*30; % RES OFFSET - SET TO 0 FOR MAGNUS TERM EFFECT STUDIES
overRotation = 0; %0 is a perfect pulse, +/- 0.01 is a 1% over/under rotation, etc.
phaseTrans = 0; %0 is no phase transient, 1 is a pi/2 phase transient

%% Configuration

% Derived parameters
dim = 2^N;
f1 = 1/4/pulse; %Can adjust f1 and w1 by changing 'pulse' variable WHAT IS F1?
w1 = 2*pi*f1;
rotationError = 1+overRotation; %1 is a perfect pulse

%Configure Test Variable
testValueMax = getVarMax(testVarName,varMaxMod); %GIVES THE MAX VALUE FOR THE TEST VAR
% VALUES HARD CODED IN THIS FUNCTION 

initVars 
initCollectiveObs

%% Initialize Hamiltonians
% 
%  Hdips = cell(couplingsCount,1);
%  % Generate [couplingsCount] dipole Hamiltonians, each with a different coupling matrix
%  for j=1:couplingsCount
%      dip = abs(randn(N));
%      dip = triu(dip,1) + triu(dip,1)';
%      Hdips{j} = getHdip(N, dim, x, y, z, dip);
%  end
load('Standard_Dipole_Hamiltonians(4).mat','Hdips'); %CHOOSE RAND OR LOAD THIS??

%% Initialize Result Storage Arrays

testVars = zeros(testValueCount,1);

% Size of Magnus Terms
results_hsizes = zeros(length(sequence_list),length(testVars),maxTerm+1);
raw_hsizes = zeros(length(sequence_list),length(testVars),couplingsCount,maxTerm+1);
traces_normalized =  zeros(length(sequence_list),length(testVars),couplingsCount,maxTerm+1);
traces_normalized_avg = zeros(length(sequence_list),length(testVars),maxTerm+1);

% % Hn Fidelities
% results_f = zeros(length(testVars),maxTerm+1); % GEN - EXP UNITARY WITH ERRORS
% results_Df = zeros(length(testVars),maxTerm+1); % D IS FOR DELTA PULSES
% raw_f = zeros(length(testVars),couplingsCount,maxTerm+1);
% raw_Df = zeros(length(testVars),couplingsCount,maxTerm+1);
% 
% % Time-Suspension Fidelities
% results_fTS = zeros(length(testVars),1); % GEN - EXP UNITARY WITH ERRORS
% results_DfTS = zeros(length(testVars),1); % D IS FOR DELTA PULSES
% raw_fTS = zeros(length(testVars),couplingsCount);
% raw_DfTS = zeros(length(testVars),couplingsCount);
% 
% % Commutation with Other Terms
% results_C0 = zeros(length(testVars),maxTerm+1); %??
% results_CS = zeros(length(testVars),maxTerm+1); %??
% raw_C = zeros(length(testVars),couplingsCount,maxTerm+1);
% raw_CS = zeros(length(testVars),couplingsCount,maxTerm+1);

%% Iterate over different parameter values, compute AHT terms and fidelities

for d=1:length(testVars) % d IS testVars INDEX
    
    % Configure Test Variable 
    if strcmp(testVarName,'tau')||strcmp(testVarName,'Tau')
        tau = pulse + d*(testValueMax/testValueCount);
        %tau = 10^(-d);
        testVars(d) = tau;
    elseif strcmp(testVarName,'Delta')||strcmp(testVarName,'delta')||strcmp(testVarName,'delta_lo')
        Delta = 2*(d-(testValueCount/2))*(testValueMax/testValueCount);
        testVars(d) = Delta;
    elseif strcmp(testVarName,'Coupling')||strcmp(testVarName,'coupling')||strcmp(testVarName,'coupling_lo')
        coupling = d*(testValueMax/testValueCount);
        testVars(d) = coupling;
    elseif strcmp(testVarName,'Overrotations')||strcmp(testVarName,'overrotations')||strcmp(testVarName,'overrot_hi')
        overRotation = 2*(d-(testValueCount/2))*(testValueMax/testValueCount);
        rotationError = 1+overRotation;
        testVars(d) = overRotation;
    elseif strcmp(testVarName,'phaseTrans')
        phaseTrans = 2*(d-(testValueCount/2))*(testValueMax/testValueCount);
        testVars(d) = phaseTrans;
    end
    
    for seq_ind = 1:length(sequence_list)
        sequenceName = sequence_list{seq_ind} ;  % select sequence to test over.
    
        % Configure Test Sequence
        sequence = getSequence(sequenceName,X,Y,Z); %X AND Y DEFINED IN initVars    
    
        % Add overrotations and phase transients
        for editPulse = 1:length(sequence.Pulses)
            sequence.Pulses{editPulse}=pulseError(sequence.Pulses{editPulse},rotationError,phaseTrans,X,Y);
        end
        
        Pulses = sequence.Pulses;%THESE ARE NOW SET TO THE DELTA PULSES WITH ERRORS??
        Taus = tau * sequence.Taus;        
        tCyc = sum(Taus);
        
        % For each dipolar coupling matrix
        for c=1:couplingsCount  % c IS couplingsCount INDEX
            
            Hdip = Hdips{c};
            %Hsys = Hdip*coupling + Z*Delta; 
            Hsys =  Z*Delta; 
            %hsys = matOrder(Hsys); 
            hsys = matOrder(Hdip*coupling); 
            %Hsys = Hdip*coupling + Z*Delta;    

            toggledHsys = {};
            for p = 0:length(Pulses)
                toggledHsys{p+1} = getURF(p)'*Hsys*getURF(p); %#ok<*SAGROW> %getURF IS FOR CALCULATING HINT IN TOGGLING FRAME
            end
            
            %Calculate Magnus terms
            knownOmegas = {};
            knownPs = {};
            knownQs = {};        
            MagnusTerms = {};
            getMagnusTerms;
            
        %     % obtain experimental unitary
        %     testUnitaries = getTestUnitaries(sequence,Hsys,pulse,tau,f1);
        %     %pulse HERE IS THE PULSE LENGTH
        %     expUnitary = testUnitaries{1};
        %     deltaUnitary = testUnitaries{2};
        % 
        %     % obtain theoretical unitaries from AHT and compute Fidelity
        %     AHTUnitaries = {};
        %     sumAH = sparse(dim,dim);
        % 
        %     for au=1:maxTerm+1
        %         sumAH = sumAH + MagnusTerms{au};
        %         AHTUnitaries{au} = expm(-1i*sumAH*2*pi*tCyc);
        %         raw_f(d,c,au)=metric(expUnitary, AHTUnitaries{au}, N);
        %         raw_Df(d,c,au) = metric(deltaUnitary, AHTUnitaries{au}, N);
        %     end
        % 
        %     % Time-suspension fidelities
        %     raw_fTS(d,c) = metric(expUnitary,speye(dim,dim),N);
        %     raw_DfTS(d,c) = metric(deltaUnitary,speye(dim,dim),N); 
        % 
         end
        
        % Average Raw Results   
        for mt = 1:maxTerm+1
            results_hsizes(seq_ind,d,mt)=mean(raw_hsizes(seq_ind,d,:,mt));
            traces_normalized_avg(seq_ind,d,mt) = mean(traces_normalized(seq_ind,d,:,mt));
            % results_f(d,mt)=mean(raw_f(d,:,mt));
            % results_Df(d,mt)=mean(raw_Df(d,:,mt)); 
        end
    
        % results_fTS(d)=mean(raw_fTS(d,:));
        % results_DfTS(d)=mean(raw_DfTS(d,:));
    
        % progress tracker for my impatient self
        strcat(testVarName,'_',sequenceName,'_',string(d))
    end
end

results_norm_hsizes_onlyDelta = results_hsizes;
trace_norm_avg_onlyDelta = traces_normalized_avg;

save(['2024-12-16_onlyDelta_Magnusterms_Iterative_Coupling_420Hz_Delta_30Hz_4spins_' ...
    '4couplings_terms0-4_20_TAUvals.mat'],'results_norm_hsizes_onlyDelta',...
    'trace_norm_avg_onlyDelta','coupling','Delta',...
    'testVars','testVarName','sequence_list','N')

% %% Make Plotting More Convenient by saving colors now
% myColors = {[230 25 75],[245 130 48], [210 245 60], [60 180 75],  [0 130 200], [145 30 180]};
% % colors are:   red          orange      lime          green           blue         violet
% myColorsPlus = {[250 190 212],[230 25 75],[245 130 48],  [210 245 60], [60 180 75], [70 240 240], [0 130 200], [145 30 180], [0 0 0]};
% % colors are:        pink          red      orange           lime          green        cyan         blue         violet       black  
% for i=1:length(myColorsPlus)
%     if i<=length(myColors)
%         myColors{i} = myColors{i}/255;
%     end
%     myColorsPlus{i} = myColorsPlus{i}/255;    
% end

% %% Save Result Output
% its = num2str(tau*10^7);
% tauString = strcat(its(1),',',its(2));
% if strcmp(testVarName,'coupling')||strcmp(testVarName,'Coupling')||strcmp(testVarName,'coupling_lo')
%     fileDescriptor = strcat(date,'_',sequenceName,'_',testVarName,tauString,']',string(Delta),'_REC_',string(maxTerm),'.mat');
% elseif strcmp(testVarName,'tau')||strcmp(testVarName,'Tau')
%     fileDescriptor = strcat(date,'_',sequenceName,'_',testVarName,string(coupling),']',string(Delta),'_REC_',string(maxTerm),'.mat');
% end
% 
% save(fileDescriptor)