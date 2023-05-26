%% Contributions by Nikhil Challa
% inputs
% Array : Contains all subcarriers stacked across consecutive time
% The function does not perform frequency correction
function [PhaseResultOriginal,PhaseResultUnwind,CorrPhase] = F_UnwrapPhaseMonotone(Array,Pos,midPoint)

    PhaseResultOriginal = Array;
    PhaseResultUnwind = Array;
    CorrPhase = zeros(1,length(Array));

    for indx = 2:length(Array)

        if indx < midPoint - 2
            if Array(indx) > Array(indx - 1)
                CorrPhase(indx) = CorrPhase(indx - 1) - 2*pi;
            else
                CorrPhase(indx) = CorrPhase(indx - 1);
            end
        elseif indx > midPoint + 2 
            if Array(indx) < Array(indx - 1)
                CorrPhase(indx) = CorrPhase(indx - 1) + 2*pi;
            else
                CorrPhase(indx) = CorrPhase(indx - 1);
            end
        else
            if Array(indx) - Array(indx - 1) > pi
                CorrPhase(indx) = CorrPhase(indx - 1) - 2*pi;
            elseif Array(indx) - Array(indx - 1) < -pi
                CorrPhase(indx) = CorrPhase(indx - 1) + 2*pi;
            else
                CorrPhase(indx) = CorrPhase(indx - 1);
            end

        end

        PhaseResultUnwind(indx) = PhaseResultUnwind(indx) + CorrPhase(indx);

    end

end