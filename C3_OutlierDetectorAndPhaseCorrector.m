%% Contributions by Nikhil Challa
clear;clc; close all;

tic 

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
    
figure(31)
tiledlayout(2,2,'Padding','Compact');
sgtitle('Phase value for each channel with phase unwrapping (Error in Channel 3)','FontSize',26) 
for ChNo = ChNoRange
    nexttile
    hold on
    plot(MeasChIndxRange,PhaseCorrectedY(ChNo+1,1:CIR_RATE:end),MeasChIndxRange,PhaseCorrectedUnwrapY(ChNo+1,1:CIR_RATE:end))
    xline(BestLogIndexes)
    yline([-pi,pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi])
    ylabel("Phase in radians")
    xlabel("Log Sample Index (every sec/skip rest)")
    set(gca,"FontSize",18)
    titleStr = "Channel No:" + string(ChNo);
    title(titleStr,'FontSize',22)
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
BadIndx = 0;
BadCh = -100;
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

figure(32)
tiledlayout(2,2,'Padding','Compact');
sgtitle('Phase value for each channel with phase unwrapping and error correction','FontSize',26) 
for ChNo = ChNoRange
    nexttile
    hold on
    plot(MeasChIndxRange,PhaseCorrectedY(ChNo+1,1:CIR_RATE:end),MeasChIndxRange,PhaseCorrUnWrapOutlierFixedY(ChNo+1,1:CIR_RATE:end))
    xline(BestLogIndexes)
    yline([-pi,pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi])
    ylabel("Phase in radians")
    xlabel("Log Sample Index (every sec/skip rest)")
    set(gca,"FontSize",18)
    titleStr = "Channel No:" + string(ChNo);
    title(titleStr,'FontSize',22)
    hold off
end

figure(33)
tiledlayout(2,2,'Padding','Compact');
sgtitle('Incorrect Phase unwrapping example','FontSize',26)
for ChNo = ChNoRange
    % We attempt to plot one (and only one) case of incorrect phase wrapping
    IndxRange = BadIndx - 200 : 1 : BadIndx;  % since we start from bad case and go backward, the issue will be seen in previous 200 samples at max
    Offset = PhaseCorrectedY(:,IndxRange(1)) - PhaseCorrectedUnwrapY(:,IndxRange(1));
    TempPh = PhaseCorrectedY(:,IndxRange) - Offset;

    % to determine limits on Y-axis
    XlimMin = floor(min([ PhaseCorrectedUnwrapY(:,IndxRange) TempPh ],[],'all')) - 1;
    XlimMax = ceil(max([PhaseCorrectedUnwrapY(:,IndxRange) TempPh ],[],'all')) + 1;

    nexttile
    hold on
    plot( IndxRange, PhaseCorrectedUnwrapY(ChNo+1,IndxRange), "rx",IndxRange, TempPh(ChNo+1,:),'bo')
    ylabel("Phase in radians")
    xlabel("Log Sample Index (sample time)")
    ylim([XlimMin  XlimMax])
    xlim([min(IndxRange) max(IndxRange)])
    legend("Unwrapped Attempt","Original Phase","Location","best")
    set(gca,"FontSize",18)
    titleStr = "Channel No:" + string(ChNo);
    title(titleStr,'FontSize',22)
    hold off

end



save('MatFiles/PhaseCorrUnWrapOutlierFixedY.mat','PhaseCorrUnWrapOutlierFixedY');
save('MatFiles/CorrPhaseOutlierFixed.mat','CorrPhaseOutlierFixed');


% Attempt to box out the incorrect phase unwrapping
% the factor of 20, needs to be decided based on number of 
window = floor(length(MeasChIndxRange)/length(BestLogIndexes));
Range = (BadIndx-CIR_RATE*window):(BadIndx+CIR_RATE*window);
Lx = floor(BadIndx/200) - window;
Ly = floor(min(PhaseCorrectedUnwrapY(BadCh,Range)));
Hx = window*2;
Hy = ceil(max(PhaseCorrectedUnwrapY(BadCh,Range))) - Ly;

figure(34)
tiledlayout(1,2,'Padding','Compact');
sgtitle('Before and After phase unwrapping error for Channel3','FontSize',26)

nexttile
hold on
plot(MeasChIndxRange,PhaseCorrectedY(BadCh,1:CIR_RATE:end),MeasChIndxRange,PhaseCorrectedUnwrapY(BadCh,1:CIR_RATE:end))
xline(BestLogIndexes)
yline([-pi,pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi])
ylabel("Phase in radians")
xlabel("Log Sample Index (every sec/skip rest)")
rectangle('Position',[Lx Ly Hx Hy],'EdgeColor','g')
set(gca,"FontSize",18)
titleStr = "Channel No:" + string(ChNo) + "  Phase value with Unwrapping Error";
title(titleStr,'FontSize',22)
hold off 

nexttile
hold on
plot(MeasChIndxRange,PhaseCorrectedY(BadCh,1:CIR_RATE:end),MeasChIndxRange,PhaseCorrUnWrapOutlierFixedY(BadCh,1:CIR_RATE:end))
xline(BestLogIndexes)
yline([-pi,pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi])
ylabel("Phase in radians")
xlabel("Log Sample Index (every sec/skip rest)")
rectangle('Position',[Lx Ly Hx Hy],'EdgeColor','g')
set(gca,"FontSize",18)
titleStr = "Channel No:" + string(ChNo) + "  Phase value with Unwrapping Error Fixed";
title(titleStr,'FontSize',22)
hold off 


toc

    








