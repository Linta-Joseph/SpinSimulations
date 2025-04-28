% Function that defines and evolves a random state and calculate the 
% entanglement entropy for a bipartite cut, given certain input parameters.

function [EE_avg] = ESS_and_EE(M,N,v,g,tstart,tinc,tfin,reals,sym,tnum,J,h)

    % Pauli matrices
    z=sparse([1 0; 0 -1]);x=sparse( [ 0 1;1 0]); y=1i*sparse([0 -1;1 0]);
    Sz=sparse(0.5.*[1 0; 0 -1]); Sx=sparse(0.5.*[ 0 1;1 0]); Sy=1i*sparse(0.5.*[0 -1;1 0]);
    ep=sparse([1 0; 0 0]);
    em=sparse([0 0; 0 1]);
    id=speye(2);p=sparse([0 1;0 0]);m=sparse([0 0; 1 0]);
    dim = 2^(M+N);
    
    % set u with or without a normalization condition
    %u = 0.25/J;
    u = sqrt(0.25-(J*v)^2)*(1/J); 
    
    assignin('base', 'M', M);
    assignin('base', 'N', N);
    assignin('base', 'dim', dim);
    assignin('base', 'dimN', 2^N);
    assignin('base', 'g', g);
    assignin('base', 'u', u);
    assignin('base', 'v', v);
    assignin('base', 'tstart', tstart);
    assignin('base', 'tinc', tinc);
    assignin('base', 'tfin', tfin);
    assignin('base', 'tnum', tnum);
    assignin('base', 'reals', reals);
    assignin('base', 'sym', sym);
    
    % Nearest neighbor couplings
    a=diag(ones(N+M-1,1),1)+diag(ones(N+M-1,1),-1);
    
    % Random Psi
    up = [1,0];
    down = [0,1];
    
    % H is the Hamiltonian that the spins evolve under
    H=sparse(dim,dim);
    
    % build up Hamiltonian
    for k=1:(N+M) 
 	    for l=k+1:(N+M)
    
	        H=H+a(k,l)*J*(((u+v)/2)*mykron(speye(2^(k-1)),Sx,speye(2^(l-k-1)),Sx,speye(2^(N+M-l)))+...
	            ((v-u)/2)*mykron(speye(2^(k-1)),Sy,speye(2^(l-k-1)),Sy,speye(2^(N+M-l)))-...
	            v*mykron(speye(2^(k-1)),Sz,speye(2^(l-k-1)),Sz,speye(2^(N+M-l)))) ;
    
	    end
    end
    
    % matrix to store EE 
    EE = zeros(reals,tnum);
     
     delete(gcp('nocreate'))
     poolobj=parpool('local',16); % start at parallel pool of 16 cores on the local node
     parfor i = 1:reals
        % temporary cell to avoid parfor issues
        eig_tmp=cell(1,tnum);
        ee_tmp=zeros(1,tnum);
        
        Htot = sparse(dim,dim);
    
        Psi_rand = 1;
        
        % product state - product of single spin product spins chosen
        % randomly from the Bloch sphere
        for j=1:(M+N)  
            a= rand;
            theta=acos(2*a-1); 
            phi=rand*2*pi; 
    
            new_rand = cos(theta).*up+exp(1i*phi)*sin(theta).*down;
            Psi_rand = mykron(Psi_rand,new_rand);
        end
    
        X = sparse(dim,dim);
        Y = sparse(dim,dim);
        Z = sparse(dim,dim);    % collective nuclear spin polarization
    
        % build up disorder term
        for k=1:(N+M)
            
            gxi = h*g*(-1+2*rand);
            X = X + gxi*mykron(speye(2^(k-1)),Sx,speye(2^(N+M-k)));
            
            % gyi = 0;%g*(-1+2*rand);
            Y = 0;%Y + gyi*mykron(speye(2^(k-1)),Sy,speye(2^(N+M-k)));
            
            gzi = h*g*(-1+2*rand);
            Z = Z + gzi*mykron(speye(2^(k-1)),Sz,speye(2^(N+M-k)));
    
        end
    
        % total Hamiltonian
        if sym == 'z'
            Htot = H + Z;
        elseif sym == 'xz'
            Htot = H + X + Z;
        elseif sym == 'xyz'
            Htot = H + X + Y + Z;
        end
        
        
        % define Unitary operators for first and subsequent time steps
        % evolve from zero or to initial time if not from zero
        if tstart == 0
            
            % define unitary for time evolution
            U=sparse(expm(-1i*Htot*tinc));
            % initial state is a product state
            Psi_t=Psi_rand';
          
        else
            
            % define unitaries for initial and subsequent time evolutions
            U0=sparse(expm(-1i*Htot*tstart));
            
            % define unitary for subsequent time steps
            U=sparse(expm(-1i*Htot*tinc));    
            
            % evolve to initial time
            Psi_t=U0*Psi_rand';
            
        end
        
        % calculate eigenvalues of reduced density matrix
        rho = Psi_t*Psi_t';
        rhoA = PartialTrace(rho,1,[2^M,2^N]);
        eig_tmp{1,1} = flip(sort(eig(rhoA)));
        
        % sample data at subsequent times
        for t = 2:tnum
            
            Psi_t=U*Psi_t;
            
            % reduced density matrix
            rho = Psi_t*Psi_t';
            rhoA = PartialTrace(rho,1,[2^M,2^N]);
            eig_tmp{1,t} = flip(sort(eig(rhoA)));
    
            %Calculate EE
            a=eig_tmp{1,t};
            b=0;
            for n = 1:length(eig_tmp{1,t})
                b = b - a(n)*log(a(n));
            end
            ee_tmp(1,t)=b;
        end
        EE(i,:) = ee_tmp;
    end
    EE_avg = mean(EE,1);
    delete(poolobj);
    clear new_rand angles theta phi1 phi2 i j k t;
end
