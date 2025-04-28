% Modified and expanded from code written by Kent Ueno
% Linta Joseph; February 2025 

% Function that defines a Hamiltonian given input parameters and performs: 
% (1) eigenvalue/vector decomposition, stores the eigenvalues; 
% (2) Schmidt decompositions for a range of eigenvectors, 
% then stores each Schmidt spectrum (zeta) as an array in a column 
% of the cell zetcell;
% (3) calculates the average entnglement entropy from 
% the schmidt decomposition for the same set of mid-spectrum eigensgtates.

function [EE_avg,eigvals,zetcell] = Schmidt_ES_NMR(M,N,u,v,ximax,yimax,zimax,J,h)

    % Pauli/spin matrices
    z=sparse([1 0; 0 -1]); x=sparse( [ 0 1;1 0]); y=1i*sparse([0 -1;1 0]);
    Sz=sparse(0.5.*[1 0; 0 -1]); Sx=sparse(0.5.*[ 0 1;1 0]); Sy=1i*sparse(0.5.*[0 -1;1 0]);
    ep=sparse([1 0; 0 0]);
    em=sparse([0 0; 0 1]);
    id=speye(2);p=sparse([0 1;0 0]);m=sparse([0 0; 1 0]);
    dim = 2^(M+N);
    dimN = 2^N;
    
    spectrum_third = round(dim/3); % a third of the energy spectrum size
    
    % assign some variables to the workspace
    assignin('base', 'M', M);
    assignin('base', 'N', N);
    assignin('base', 'dim', dim);
    assignin('base', 'dimN', dimN);
    assignin('base', 'u', u);
    assignin('base', 'v', v);
    assignin('base', 'yimax', yimax);
    assignin('base', 'zimax', zimax);
    assignin('base', 'zimax', zimax);
    assignin('base', 'J', J);
    assignin('base', 'h', h);
    
    % H is the Hamiltonian that the spins evolve under
    H=sparse(dim,dim);
    
    % nearest neighbor couplings
    a=sparse(diag(ones(N+M-1,1),1)+diag(ones(N+M-1,1),-1));
    
    % build up spin chain properties - 10-28-2024 -- REPLACING PAULI OPS WITH
    % SPIN OPS
    for k=1:(N+M) 
 	    for l=k+1:(N+M)
    
	        H=H+a(k,l)*J*(((u+v)/2)*mykron(speye(2^(k-1)),Sx,speye(2^(l-k-1)),Sx,speye(2^(N+M-l)))+...
	            ((v-u)/2)*mykron(speye(2^(k-1)),Sy,speye(2^(l-k-1)),Sy,speye(2^(N+M-l)))-...
	            v*mykron(speye(2^(k-1)),Sz,speye(2^(l-k-1)),Sz,speye(2^(N+M-l)))) ;
    
	    end
    end
    
    % collective nuclear spin polarization
    X = sparse(dim,dim);
    Y = sparse(dim,dim);
    Z = sparse(dim,dim);
    
    % build up disorder terms - 10-28-2024 -- REPLACING PAULI OPS WITH
    % SPIN OPS
    for k=1:(N+M)
    
        xi = h*ximax*(-1+2.0*rand); %random x field strength 
        % build up x disorder
	    X = X + xi.*mykron(speye(2^(k-1)),Sx,speye(2^(N+M-k)));
    
	    % yi = yimax*(-1+2.0*rand);
        yi = 0;%h*yimax*(-1+2.0*rand);
	    Y = 0;%Y + yi.*mykron(speye(2^(k-1)),Sy,speye(2^(N+M-k)));
    
	    % zi = zimax*(-1+2.0*rand);
        zi = h*zimax*(-1+2.0*rand);  %random z field strength
        % build up z disorder
        Z = Z + zi.*mykron(speye(2^(k-1)),Sz,speye(2^(N+M-k)));
    
    end
    
    % total Hamiltonian
    H = sparse(H + X + Y + Z);
    
    %check Hermitian 
    % ishermitian(H)
    
    % sort eigenenergies and eigenvectors in ascending order
    [V,D] = eig(full(H));
    [~,I] = sort(diag(D));
    V = V(:,I);
    
    eigvals{1,1} = sort(diag(D));
    
    clear D I;
    
    % middle third of eigenenergy spectrum
    V = V(:,spectrum_third:(2*spectrum_third));
    
    % matrix to store each array of zeta values
    zetcell = cell(1,min(size(V)));
    EE = zeros(1,min(size(V)));
    
    % Schmidt decomposition of eigenvectors
    for c=1:(min(size(V)))
        
	    S = single(SchmidtDecomposition(V(:,c),[2^(M),2^(N)]));
        % Calculate entanglement entropy
        EE(1,c)=-sum((S.^2).*log(S.^2));
        % Store each zeta array as a column of the zeta cell
	    zetcell{1,c} = (-2).*log(S);
    end
    % Average of the entanglement entropy over the mid-spectrum eigenstates
    EE_avg = mean(EE,2);
    clear V S;

end