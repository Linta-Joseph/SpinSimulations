% Arguments: %1) sequenceName in a specific capital letter acronym format
% 2) X and Y spin operators. 
% Returns Struct sequence with arrays of Pulses and Taus, 

function sequence = getSequence(sequenceName, X, Y)

    %WAHUHA
    if strcmp(sequenceName, 'WHH')
        sequence.Pulses = {-X, Y, -Y, X}; %(rev ppr, br ppr, type 1b)
        sequence.Taus = [1 2 1 2];

    %MREV-8
    elseif strcmp(sequenceName, 'MREV8')
        sequence.Pulses = {X, Y, -Y, -X, -X, Y, -Y, X}; % (rev ppr version 1)
        sequence.Taus = [1 2 1 2 1 2 1 2];

    %MREV-16
    elseif strcmp(sequenceName, 'MREV16')
        sequence.Pulses = {-X, Y, -Y, X, X, Y, -Y, -X, -X, -Y, Y, X, X, -Y, Y, -X}; % (Ladd ppr, without pi pulses)
        sequence.Taus = [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2];
        
    %BR-24    
    elseif strcmp(sequenceName, 'BR24')
        sequence.Pulses = {X, Y, -Y, -X, -X, Y, -Y, X, Y, X, -X, -Y, -Y, X,...
            Y, X,-X, -Y, -Y, X,-X, Y, -X, Y }; % Burum, Rhim, analysis 3 (1a1b2a(-yx)2a2b(-xy))
        sequence.Taus = [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2];
        
    %CORY 48
    elseif strcmp(sequenceName, 'CORY48')
        sequence.Pulses = {X, Y, -X, Y, X, Y, X, Y, X, -Y, X, Y, -Y, -X, Y, ...
            -X, -Y, -X, -Y, -X, -Y, X, -Y, -X, -X, Y, -X, -Y, -X, Y, X, -Y, ...
            -X, -Y, X, -Y, Y, -X, Y, X, Y, -X, -Y, X, Y, X, -Y, X}; % cory 1990
        sequence.Taus = [1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 ... 
            2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2];
        
    %YXX-48
    elseif strcmp(sequenceName,'YXX48')
        sequence.Pulses = {Y, -X,-X,Y,-X,-X,-Y,X,X,Y,-X,-X,-Y,X,X,-Y,X,X,Y, ...
            -X,-X,Y,-X,-X,-Y,X,X,Y,-X,-X,-Y,X,X,-Y,X,X,Y,-X,-X,-Y,X,X,Y,-X, ...
            -X,-Y,X,X}; %Peng RL ppr
        sequence.Taus = ones(48,1);

    
    %YXX-24
    elseif strcmp(sequenceName,'YXX24')
        sequence.Pulses = {-Y,X,-X,Y,-X,-X,Y,-X,X,-Y,X,X,Y,-X,X,-Y,X,X,-Y, ... 
            X,-X,Y,-X,-X}; %Peng RL ppr
        sequence.Taus = ones(24,1);

       
    %SED-48
    elseif strcmp(sequenceName,'SED48')
        sequence.Pulses = {-X, -Y, -X, Y, -X, Y, -Y, -X, Y, -X, -Y, -X, X, ... 
            Y, X, -Y, X, -Y, Y, X, -Y, X, Y, X, X, Y, Y, X, Y, X, -Y, X, -X, ...
            Y, -X, -Y, -X, -Y, -X, Y, -Y, X, -Y, -X, -Y, -X, X, Y}; %Owen thesis    
        sequence.Taus=[1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, ...
            1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, ...
            2, 1, 2, 1, 2, 1, 2, 1, 2];

    %SED-36
    elseif strcmp(sequenceName,'SED36')
        sequence.Pulses ={X, Y, X, -Y, Y, -X, -Y,-X, -Y, -X, X, Y, Y, -X, Y, ...
            -X, Y, X, -X, -Y, X, -Y, X, -Y};
        sequence.Taus=[1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, ...
            1, 2, 1, 2, 1, 2]; % Madhu e-mail (01/05/24)
    
    %SED-36err
    elseif strcmp(sequenceName,'SED36err')
        sequence.Pulses =[-Y, X, -Y, -X, -Y, X, -X, -Y, -X, Y, X, Y, -Y, X, ... 
            Y, X, -Y, X, -X, -Y, X, -Y, X, Y];
        sequence.Taus=[1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, ...
            1, 2, 1, 2, 1, 2]; % Madhu e-mail (01/05/24)
        
    end
end
