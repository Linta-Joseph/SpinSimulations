% GENERATES N-SPIN DIPOLAR HAMILTONIAN Hdip 
% COUPLINGS SPECIFIED BY THE SYMMETRIC MATRIX a

function Hdip = getHdip(N, dim, x, y, z, a)

% N: number of spins
% x/y/z: Spin matrices (0.5*sigma)
% a: dipolar interaction strengths (NxN symmetric matrix)

Hdip = sparse(dim, dim);

for k=1:N
    for h=k+1:N
        % Hdip = a_{kh} (3Z - X - Y)
        Hdip=Hdip+a(k,h)*(2*mykron(speye(2^(k-1)),z,speye(2^(h-k-1)),z,speye(2^(N-h)))-... 
            mykron(speye(2^(k-1)),x,speye(2^(h-k-1)),x,speye(2^(N-h)))-...
            mykron(speye(2^(k-1)),y,speye(2^(h-k-1)),y,speye(2^(N-h)))) ;             
    end
end

end

%MULTIPLIED BY Omega_D LATER. BUT, ARE THERE OTHER MISSING FACTORS?