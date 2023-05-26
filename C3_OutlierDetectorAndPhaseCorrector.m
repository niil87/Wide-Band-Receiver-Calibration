%% Contributions by Nikhil Challa
clear;clc; close all;

% The code has the following limitations
% 1. It assumes only 1 outlier at a time for input signal in array of antennas
% 2. It assumes 3 or more array of antennas otherwise outlier will not work
% Note that this is very rare case and should not be encounter if the
% antenna is moved slowly hence the limitation.

BestLogIndexes = load('MatFiles/BestLogIndexes.mat').BestLogIndexes;
PhaseCorrectedY = load('MatFiles/PhaseCorrectedY.mat').PhaseCorrectedY;
PhaseCorrectedUnwrapY = load('MatFiles/PhaseCorrectedUnwrapY.mat').PhaseCorrectedUnwrapY;
CorrPhase = load('MatFiles/CorrPhase.mat').CorrPhase;
CIR_RATE = 200;
ChNoRange = 0:3;
MeasChIndxRange = 11:635; 
PhaseCorrUnWrapOutlierFixedY = PhaseCorrectedUnwrapY;
CorrPhaseOutlierFixed = CorrPhase;
for ChNo = ChNoRange

    figure(1 + ChNo)
    hold on
    plot(MeasChIndxRange,PhaseCorrectedY(ChNo+1,1:CIR_RATE:end),MeasChIndxRange,PhaseCorrectedUnwrapY(ChNo+1,1:CIR_RATE:end))
    xline(BestLogIndexes)
    yline([-pi,pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi])
    ylabel("Phase in radians")
    xlabel("Samples/CIR-RATE")
    titleStr = "Channel No:" + string(ChNo) + "  Phase value after removing frequency offset";
    title(titleStr)
    hold off

end


% Ideally we would like to check every index but the jitters can cause a
% lot of confusion on where the false phase jump happened. Instead if we
% check every 200 indexes (similar one use for plotting) we will not see
% these jitters and the accuracy on false phase jump detection is lot
% higher

outlier = zeros(1,floor(length(PhaseCorrUnWrapOutlierFixedY)/200));

outlier(1) = ~ isempty( find(isoutlier(PhaseCorrUnWrapOutlierFixedY(:,1)),1) );
outlier(2) = ~ isempty( find(isoutlier(PhaseCorrUnWrapOutlierFixedY(:,201)),1) );
outIndx = 3;
for indx = 401 : 200 : length(PhaseCorrUnWrapOutlierFixedY)
    outlier(outIndx) = ~ isempty( find(isoutlier(PhaseCorrUnWrapOutlierFixedY(:,indx)),1) );
    if (outlier(outIndx - 2) && outlier(outIndx - 1) && outlier(outIndx))
        BadIndx = (outIndx - 3)*200+1;
        BadCh = find(isoutlier(PhaseCorrUnWrapOutlierFixedY(:,BadIndx)),1);
        GoodCh = setdiff(1:4,BadCh);
        BaseCh = GoodCh(1);
        PhDiff = PhaseCorrUnWrapOutlierFixedY(BaseCh,BadIndx) - PhaseCorrUnWrapOutlierFixedY(BadCh,BadIndx);
        
        % Note that phase unwrap are always multiples of 2pi
        PhShift = 0;
        if PhDiff > 0
            PhShift = 2*pi;
        else 
            PhShift = -2*pi;
        end
        disp(["Outlier Detected at indx:" + BadIndx + " on Channel:" + string(BadCh-1) + " with Ph diff of:" + PhShift])
        
        for phIndx = BadIndx:length(PhaseCorrUnWrapOutlierFixedY)
            PhaseCorrUnWrapOutlierFixedY(BadCh,phIndx) = PhaseCorrUnWrapOutlierFixedY(BadCh,phIndx) + PhShift;
            CorrPhaseOutlierFixed(BadCh,phIndx) = CorrPhaseOutlierFixed(BadCh,phIndx) + PhShift;
        end

        % readjusting indx back to recompute outliers
        indx = BadIndx;
        outIndx = outIndx - 3;
    end
    outIndx = outIndx + 1;

end


for ChNo = ChNoRange

    figure(5 + ChNo)
    hold on
    plot(MeasChIndxRange,PhaseCorrectedY(ChNo+1,1:CIR_RATE:end),MeasChIndxRange,PhaseCorrUnWrapOutlierFixedY(ChNo+1,1:CIR_RATE:end))
    xline(BestLogIndexes)
    yline([-pi,pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi])
    ylabel("Phase in radians")
    xlabel("Samples/CIR-RATE")
    titleStr = "Channel No:" + string(ChNo) + "  Phase value after removing frequency offset and Unwrap correction";
    title(titleStr)
    hold off

end

save('MatFiles/PhaseCorrUnWrapOutlierFixedY.mat','PhaseCorrUnWrapOutlierFixedY');
save('MatFiles/CorrPhaseOutlierFixed.mat','CorrPhaseOutlierFixed');


    








