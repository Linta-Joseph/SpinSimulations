
%Arguments: %1)sequenceName in a specific capital letter acronym format
%2)X and Y spin operators. 
% Returns Struct sequence with 
% 1)arrays of Pulses and Taus, 
% 2)no: of cycles the sequence is to be applied (256 or 128)
% 3)initial density matrix -- (normal to eff. field for spectro, {X,Y,Z} for timesus)

%% NOTE: First/last 2tau spacing split to 1-1 at beginning and end -- 2024-10-02
function sequence = getSequence(sequenceName, X, Y, Z)
    %global X Y %#ok<GVMIS> 
    %WAHUHA
    if strcmp(sequenceName, 'WHH')
        sequence.Pulses = {-X, Y, -Y, X}; %(rev ppr, br ppr, type 1b)
        sequence.Taus = [1 1 2 1 1];
        sequence.num_cycle = 50;
        %sequence.rho_in ={(1/sqrt(6))*(X+Y-2.*Z)}; % effective field: 1/3(X+Y+Z)
        sequence.rho_in = {X Y Z};
        sequence.fitpar = 1*10^-3 ; 

    %MREV-8
    elseif strcmp(sequenceName, 'MREV8')
        sequence.Pulses = {X, Y, -Y, -X, -X, Y, -Y, X}; % (rev ppr version 1)
        sequence.Taus = [1 1 2 1 2 1 2 1 1];
        sequence.num_cycle = 50;
        %sequence.rho_in = {Y} ; % effective field: 1/3(X+Z)
        sequence.rho_in = {X Y Z};
        sequence.fitpar = 1*10^-3;

    %MREV-16
    elseif strcmp(sequenceName, 'MREV16')
        sequence.Pulses = {-X, Y, -Y, X, X, Y, -Y, -X, -X, -Y, Y, X, X, -Y, Y, -X}; % (Ladd ppr, without pi pulses)
        sequence.Taus = [1 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 1];
        sequence.num_cycle = 50;
        %sequence.rho_in = {Y} ; % effective field: 1/sqrt(2)(X+Z)
        sequence.rho_in = {X Y Z};
        sequence.fitpar = 1*10^-3;

        
    %BR-24    
    elseif strcmp(sequenceName, 'BR24')
        sequence.Pulses = {X, Y, -Y, -X, -X, Y, -Y, X, Y, X, -X, -Y, -Y, X,...
            Y, X,-X, -Y, -Y, X,-X, Y, -X, Y }; % Burum, Rhim, analysis 3 (1a1b2a(-yx)2a2b(-xy))
        sequence.Taus = [1 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 1];
        sequence.num_cycle = 50;
        %sequence.rho_in = {(1/3)*(X+Y-2.*Z)};?????
        sequence.rho_in = {X Y Z};
        sequence.fitpar = 1*10^-3;
        
    %CORY 48
    elseif strcmp(sequenceName, 'CORY48')
        sequence.Pulses = {X, Y, -X, Y, X, Y, X, Y, X, -Y, X, Y, -Y, -X, Y, ...
            -X, -Y, -X, -Y, -X, -Y, X, -Y, -X, -X, Y, -X, -Y, -X, Y, X, -Y, ...
            -X, -Y, X, -Y, Y, -X, Y, X, Y, -X, -Y, X, Y, X, -Y, X}; % cory 1990
        sequence.Taus = [1 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 ... 
            2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 1];
        sequence.num_cycle = 50;
        sequence.rho_in = {X Y Z};
        sequence.fitpar = 10*10^-3;
        
    %YXX-48
    elseif strcmp(sequenceName,'YXX48')
        sequence.Pulses = {Y, -X,-X,Y,-X,-X,-Y,X,X,Y,-X,-X,-Y,X,X,-Y,X,X,Y, ...
            -X,-X,Y,-X,-X,-Y,X,X,Y,-X,-X,-Y,X,X,-Y,X,X,Y,-X,-X,-Y,X,X,Y,-X, ...
            -X,-Y,X,X}; %Peng RL ppr
        sequence.Taus = ones(49,1);
        sequence.Taus(1)=0;
        sequence.num_cycle = 50;
        sequence.rho_in = {X Y Z};
        sequence.fitpar = 10*10^-3;
    
    %YXX-24
    elseif strcmp(sequenceName,'YXX24')
        sequence.Pulses = {-Y,X,-X,Y,-X,-X,Y,-X,X,-Y,X,X,Y,-X,X,-Y,X,X,-Y, ... 
            X,-X,Y,-X,-X}; %Peng RL ppr
        sequence.Taus = ones(25,1);
        sequence.Taus(1)=0;
        sequence.num_cycle = 50;
        sequence.rho_in = {X Y Z};
        sequence.fitpar = 10*10^-3;

       
    %SED-48
    elseif strcmp(sequenceName,'SED48')
        sequence.Pulses = {-X, -Y, -X, Y, -X, Y, -Y, -X, Y, -X, -Y, -X, X, ... 
            Y, X, -Y, X, -Y, Y, X, -Y, X, Y, X, X, Y, Y, X, Y, X, -Y, X, -X, ...
            Y, -X, -Y, -X, -Y, -X, Y, -Y, X, -Y, -X, -Y, -X, X, Y}; %Owen thesis    
        sequence.Taus=[1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, ...
            1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, ...
            2, 1, 2, 1, 2, 1, 2, 1, 1];

    %SED-36
    elseif strcmp(sequenceName,'SED36')
        sequence.Pulses ={X, Y, X, -Y, Y, -X, -Y,-X, -Y, -X, X, Y, Y, -X, Y, ...
            -X, Y, X, -X, -Y, X, -Y, X, -Y};
        sequence.Taus=[1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, ...
            1, 2, 1, 2, 1, 1]; % Madhu e-mail (01/05/24)
    
    %SED-36err
    elseif strcmp(sequenceName,'SED36err')
        sequence.Pulses ={-Y, X, -Y, -X, -Y, X, -X, -Y, -X, Y, X, Y, -Y, X, ... 
            Y, X, -Y, X, -X, -Y, X, -Y, X, Y};
        sequence.Taus=[1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, ...
            1, 2, 1, 2, 1, 1]; % Madhu e-mail (01/05/24)

    % %SED-24
    % elseif strcmp(sequenceName,'SED24')
    %     sequence.Pulses = {X, Y, -X, -Y, Y, -X, Y, -X, Y, X, -X, Y, X, Y, -X, -Y, Y, X, -Y, -X, -Y, -X, X, -Y};       
    %     sequence.Taus=[1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1];
    %     sequence.num_cycle = 50;
    %     sequence.rho_in = {X Y Z};
    %     sequence.fitpar = 10*10^-3;


    % %BR-24-pi    
    % elseif strcmp(sequenceName, 'BR24-pi')
    %     sequence.Pulses = {X, Y, -Y, -X, -X, Y, -Y, X, Y, X, -X, -Y, -Y, X, Y, X,-X, -Y, -Y, X,-X, Y, -X, Y, X*pi}; % Burum, Rhim, analysis 3 (1a1b2a(-yx)2a2b(-xy))
    %     sequence.Taus = [ 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 2 1 1 1];
    %     sequence.num_cycle = 50;
    %     %sequence.rho_in = {(1/3)*(X+Y-2.*Z)};?????
    %     sequence.rho_in = {X Y Z};
    %     sequence.fitpar = 1*10^-3;

    %DROID60
    %DIRAC
    %DROID R2-D2

        
    % %SED-96
    % elseif strcmp(sequenceName,'SED96')
    %     sequence.Pulses =  {-X, Y, -Y, -X, Y, -X, X, Y, -Y, X, -X, -Y, X, -Y, X, -Y, Y, -X, -X, -Y, Y, X, Y, X, Y, X, Y, -X, -X, Y, -Y, -X, -Y, X, X, -Y, X, Y, -X, -Y, X, -Y, Y, -X, -Y, -X, X, ...
    %         Y, Y, X, X, Y, -X, Y, X, Y, Y, X, Y, -X, -Y, -X, X, -Y, -X, -Y, X, -Y, -Y, -X, -Y, X, Y, -X, X, Y, X, Y, -X, Y, -Y, -X, -Y, X, -Y, -X, -Y, X, -X, Y, -X, -Y, -X, -Y, Y, X};   
    %     sequence.Taus=[1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, ...
    %         2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1];    
    %     sequence.num_cycle = 50;
    %     sequence.rho_in = {X Y Z};
    %     sequence.fitpar = 10*10^-3;
        
%     %AZ-48
%     elseif strcmp(sequenceName,'AZ48')
%         sequence.Pulses = {-X,Y,Y,X,Y,Y,-Y,X,X,-Y,X,X,Y,X,X,-Y,X,X,-Y,X,-Y,X,X,-Y,-X,-X,Y,Y,-X,Y,Y,-Y,X,-Y,-Y,X,-Y,X,X,-Y,X,X,-Y,-X,-X,-Y,-X,-X};
%         sequence.Taus = ones(49,1);
%         sequence.Taus(1)=0;
%         sequence.num_cycle =
%         sequence.rho_in = 
%         
%     %Symmetrized 48 from YXX-24
%     elseif strcmp(sequenceName,'YXX24S')
%         sequence.Pulses = {-Y,X,-X,Y,-X,-X,Y,-X,X,-Y,X,X,Y,-X,X,-Y,X,X,-Y,X,-X,Y,-X,-X,X,X,-Y,X,-X,Y,-X,-X,Y,-X,X,-Y,-X,-X,Y,-X,X,-Y,X,X,-Y,X,-X,Y};
%         sequence.Taus = [0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1];
%         sequence.num_cycle =
%         sequence.rho_in = 
%         
%     elseif strcmp(sequenceName,'AS10M')
%         sequence.Pulses = {-X,Y,2*Y,-Y,-X,X,Y,2*Y,-Y,X};
%         sequence.Taus = ones(11,1);
%         sequence.Taus(6) = 2;
%         sequence.num_cycle =
%         sequence.rho_in = 
%         
%     elseif strcmp(sequenceName,'AS10E')
%         sequence.Pulses = {X, -Y, Y, -X, 2*Y, X, Y, -Y, -X, 2*Y}; %(added pi-pulse with 0tau after to fix basis)
%         sequence.Taus = ones(11,1);
%         sequence.Taus(11) = 0;
%         sequence.num_cycle =
%         sequence.rho_in = 
%         
%     elseif strcmp(sequenceName,'AS10E-2cycle')
%         sequence.Pulses = {X, -Y, Y, -X, 2*Y, X, Y, -Y, -X, X, -Y, Y, -X, 2*Y, X, Y, -Y, -X};
%         sequence.Taus = ones(19,1);
%         sequence.num_cycle =
%         sequence.rho_in = 
%         
%     elseif strcmp(sequenceName,'WHHWHH')
%         sequence.Pulses = {X, -Y, Y, -X, X, -Y, Y, -X};
%         sequence.Taus = [1 1 2 1 2 1 2 1 1];
%         sequence.num_cycle =
%         sequence.rho_in = 
%             
%     elseif strcmp(sequenceName,'ML10')
%         sequence.Pulses = {X,X,Y,X,X,-X,-X,-Y,-X,-X};
%         sequence.Taus = [1 1 1 1 1 2 1 1 1 1 1];
%         sequence.num_cycle =
%         sequence.rho_in = 
%     
%     elseif strcmp(sequenceName,'MG8')
%         sequence.Pulses = {X,Y,Y,X,X,Y,Y,X};
%         sequence.Taus = [1 1 2 1 2 1 2 1 1];
%         sequence.num_cycle =
%         sequence.rho_in = 
%         
%     elseif strcmp(sequenceName,'I24')
%         sequence.Pulses = {X, Y, Y, X, -Y, -X, -Y, -X, Y, X, Y, X, -Y, -X, -X, -Y, -X, Y, -X, -Y, X, Y, -X, Y};
%         sequence.Taus = [1, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 2, 1, 1];
%         sequence.num_cycle =
%         sequence.rho_in = 
        
    end
end
