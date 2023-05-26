%% Contributions by Nikhil Challa
% inputs
% NumSC : Number of subcarriers of OFDM Signal
% Array : Contains all subcarriers stacked across consecutive time
% instances
% Option : This option is to allow you to provide phase value directly
% 1. you provide original ofdm signal
% 2. you provide phase per time stamp 
% The function does not perform frequency correction
function [PhaseResultOriginal,PhaseResultUnwind,CorrPhase] = F_UnwrapPhase(Array,NumSC,opt)

    if not(isrow(Array) || iscolumn(Array))
        disp("The signal input must be either a column or row vector and not a matrix");
        return;
    elseif isrow(Array)
        Array = transpose(Array);
    end
    NumTimeSteps = floor(length(Array)/NumSC);

    ConstPhaseShift = 0;
    PhaseResultOriginal = zeros(NumTimeSteps,1);
    PhaseResultUnwind = zeros(NumTimeSteps,1);
    CorrPhase = zeros(NumTimeSteps,1);

    if opt == 1
        Y_base = Array(1:NumSC);
        PrevAngle = 0;
    elseif opt == 2
        PrevAngle = (Array(1));
        PhaseResultOriginal(1) = PrevAngle;
        PhaseResultUnwind(1) = PrevAngle;
    end

    for timeIndx = 2:NumTimeSteps

        if opt == 1
            StrIndx = (timeIndx - 1)*NumSC + 1;
            EndIndx = (timeIndx)*NumSC;
            Y = Array(StrIndx:EndIndx);
            AngleWrtBase = angle((Y_base' *Y)/norm(Y_base)^2);
        elseif opt == 2
            AngleWrtBase = (Array(timeIndx));
        end

        PhaseResultOriginal(timeIndx) = AngleWrtBase;
        PhaseResultUnwind(timeIndx) = AngleWrtBase + ConstPhaseShift;

        if abs(AngleWrtBase - PrevAngle) > 1.5*pi && AngleWrtBase < PrevAngle
            PhaseResultUnwind(timeIndx) = PhaseResultUnwind(timeIndx) + 2*pi;
            ConstPhaseShift = ConstPhaseShift + 2*pi;
        elseif abs(AngleWrtBase - PrevAngle) > 1.5*pi && AngleWrtBase > PrevAngle
            PhaseResultUnwind(timeIndx) = PhaseResultUnwind(timeIndx) - 2*pi;
            ConstPhaseShift = ConstPhaseShift - 2*pi;
        end
        CorrPhase(timeIndx) = ConstPhaseShift;
        PrevAngle = AngleWrtBase;

    end

end