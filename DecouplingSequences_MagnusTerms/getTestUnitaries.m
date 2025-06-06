function testUnitaries = getTestUnitaries(sequence,Hsys,pulse,tau,f1)
    global dim 
    
    UtauCell = {speye(dim,dim),expm(-1i*Hsys*tau*2*pi), expm(-1i*Hsys*2*tau*2*pi)};    

    testUnitary = speye(dim,dim); % WARNING: MAY NOT FUNCTION PROPERLY AT THIS TIME \?????
    deltaUnitary = speye(dim,dim);
    
    for p=1:length(sequence.Pulses) % CHECK THE TAUS AND PULSES
        Utau = UtauCell{sequence.Taus(p)+1};
        nextUs = getNextUs(sequence,Hsys,pulse,f1,p);
        nextU = nextUs{1};
        nextUD = nextUs{2};
        testUnitary = nextU * Utau * testUnitary;
        deltaUnitary = nextUD *Utau * deltaUnitary;
    end
    testUnitary = UtauCell{sequence.Taus(length(sequence.Taus))+1} * testUnitary;
    deltaUnitary = UtauCell{sequence.Taus(length(sequence.Taus))+1} * deltaUnitary;
    
    testUnitaries = {testUnitary,deltaUnitary};
end

function nextUs = getNextUs(sequence,Hsys,pulse,f1,p)
    nextU = expm(-1i*2*pi*(Hsys+f1*sequence.Pulses{p})*pulse);% Hint EVOLUTION IS CONSIDERED DURING THE FINITE PULSE
    %nextUD
    %=expm(-1i*Hsys*2*pi*pulse/2)*expm(-1i*pi*sequence.Pulses{p}/2)*expm(-1i*Hsys*2*pi*pulse/2);
    nextUD = expm(-1i*pi*sequence.Pulses{p}/2); %DELTA PULSE. BUT CONTAINS THE ERRORS
    nextUs = {nextU,nextUD};
end

