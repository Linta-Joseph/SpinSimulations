function testUnitaries = getTestUnitaries(sequence,Hsys,pulse,tau,f1,dim,overrot,trans)%,Pulses,X)

    UtauCell = {speye(dim,dim),expm(-1i*Hsys*tau), expm(-1i*Hsys*2*tau)};    


    testUnitary = speye(dim,dim); 
    deltaUnitary = speye(dim,dim);
    
    for p=1:length(sequence.Pulses) % CHECK THE TAUS AND PULSES
        Utau = UtauCell{sequence.Taus(p)+1};
        nextUs = getNextUs(sequence,Hsys,pulse,f1,p,overrot,trans);%,Pulses,X);
        nextU = nextUs{1};
        nextUD = nextUs{2};
        testUnitary =  Utau * nextU * testUnitary; % pulse first, delay next;  no delay at the beginning of the sequence
        deltaUnitary = Utau * nextUD * deltaUnitary;
    end
    
    testUnitaries = {testUnitary,deltaUnitary};
end

function nextUs = getNextUs(sequence,Hsys,pulse,f1,p,overrot,trans)%,Pulses,X) %Pulses are original pulses wihtout errs

    nextU = expm(-1i*(pi/2)*trans*sequence.Phtr_array{p})...%trailing edge transient
            *expm(-1i*(Hsys*pulse+(pi/2)*(1+overrot)*sequence.Pulses{p}))...%finite (pi/2) pulse with common overrotation for all spins
            *expm(-1i*(pi/2)*trans*sequence.Phtr_array{p});%leading edge transient
    % Hint EVOLUTION IS CONSIDERED DURING THE FINITE PULSE

    nextUD = expm(-1i*(pi/2)*trans*sequence.Phtr_array{p})... 
            *expm(-1i*pi*(1+overrot)*sequence.Pulses{p}/2)... 
            *expm(-1i*(pi/2)*trans*sequence.Phtr_array{p}); 
    %DELTA PULSE. BUT CONTAINS THE ERRORS

    nextUs = {nextU,nextUD};
end

