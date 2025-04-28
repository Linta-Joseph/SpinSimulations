% which phase transient to associate with which pulse 
% (note that alpha can be negative)
function trP = ComputeTransients(pulseIn,X,Y)  
    if pulseIn == X
        trP = Y;
    elseif pulseIn == Y
        trP = -X;
    elseif pulseIn == -X
        trP = -Y;
    elseif pulseIn == -Y
        trP = X;
    %elseif pulseIn == X*pi %X*pi denotes pi pulse 
        %trP = roterror*X+trans*Y; %% for pi pulse; only X pulse pulses considered. 
    end
end
