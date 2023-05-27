%% Contributions by Nikhil Challa
clear;clc; close all;

tic
%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NumCenters = 20;

% this is the % of time the system was stationary (acceptance rate)
% the algo with become stricter if the sample acceptance rate is higher
% than this value
% Given you have N secs of logs, what percentage of time or logs are in
% stationary condition?
RatioStableVsMovement = 0.5;  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       

% Maybe you can provide this directly         
PhaseCorrectedUnwrapY = load('MatFiles/PhaseCorrectedUnwrapY.mat').PhaseCorrectedUnwrapY;
CIR_RATE = 200;
MeasChIndxRange = 11:635;

% We can roughly estimate the upper limit of the windowsize in secs based
% on the provided parameters assuming the person has roughly spent the same
% amount of time in each position.
UpperlimitWindowSize = floor( length(PhaseCorrectedUnwrapY) / CIR_RATE / NumCenters * RatioStableVsMovement)


% this defines how many secs of logs will we consider as 1 window
% We have selected odd numbers so its easier to pic middle log as good
% time point and code runs a lot faster
TimeSize = 1:2:UpperlimitWindowSize; 

NumCh = size(PhaseCorrectedUnwrapY,1);
LogIndexes = zeros(NumCh,NumCenters);

figure(20)
plot(PhaseCorrectedUnwrapY(1,:))
ylabel("Phase in radians")
xlabel("Log Sample Index (every sec/skip rest)")
set(gca,"FontSize",14)
title("Unwrapped phase for Channel 0",'FontSize',18)

exit = 0;
SumDLimit = 1;
while(exit == 0)
    exit = 1;

    for ChIndx = 1:NumCh
        MaxSumD = zeros(length(TimeSize),1);
        sampleIndx = 1;
        PhaseY = PhaseCorrectedUnwrapY(ChIndx,:);
        
        Xstore = [];
        Cstore = [];
        ArraySize = [];
        Xorigstore = [];

        for SizeInSec = TimeSize
            % disp(["Attempting for Window Size",SizeInSec]);
            Window = CIR_RATE*SizeInSec;
            count = 1;
            Mean = [];
            Var = [];
            while(1)
                % Going with sliding window for better accuracy
                StartIndx = CIR_RATE*(count-1) + 1;
                EndIndx = StartIndx + Window - 1;
                if EndIndx > length(PhaseY)
                    break;
                end
                
                Sample = PhaseY(StartIndx:EndIndx);
                Mean(count) = mean(Sample);
                Var(count) = var(Sample);
                count = count + 1;
            end
        
            X = zeros(count-1,3);
            X(:,1) = Mean;
            X(:,2) = Var;
            X(:,3) = 1:length(Mean);
            sumD = ones(NumCenters,1);
            Xorigstore{sampleIndx} = X;
            
            % We decrement SumDLimit by factor of 10 everytime we exceed
            % acceptance percentage limit (1 -> 0.1 -> 0.01 -> ...)
            while(max(sumD) > SumDLimit)
                [idx,C,sumD] = kmeans(X(:,1),NumCenters);
                [~,indx] = max(X(:,2));
                X(indx,:) = [];
            end
    
            if length(X(:,1)) > length(Mean)*RatioStableVsMovement
                % what it means is the number of rejected points is not
                % enough to meet the known knowledge of what % is expected to
                % be rejected! Reducing sumD requirement to make it more
                % agressive at rejection
                disp(["Acceptance percentage:" + length(X(:,1))/length(Mean) ]);
                SumDLimit = SumDLimit/10;
                disp(["Setting SumDLimit:" + SumDLimit + " and restarting computing" ]);
                exit = 0;
                break; % breaks out of SizeInSec for-loop
            end

            disp(["(Successfully below RatioStableVsMovement) Acceptance percentage:" + length(X(:,1))/length(Mean) ]);
            
            Xstore{sampleIndx} = X;
            Cstore{sampleIndx} = C;
            ArraySize(sampleIndx) = length(X(:,1));
            sampleIndx = sampleIndx + 1;
        end

        ArraySize
        
        if (exit == 0)
            break;  % breaks out of ChIndx for-loop
        end

        [~,bestIndx] = max(ArraySize);
        X = Xstore{bestIndx};
        C = Cstore{bestIndx};
        
        figure(ChIndx)
        hold on
        scatter(X(:,1),X(:,2),'.');
        scatter(C(:,1),zeros(length(C),1),'o');
        titleStr = "Channel No:" + string((ChIndx - 1)) + "  Clustering Data";
        title(titleStr);
        xlabel 'Mean'; 
        ylabel 'Variance';
        hold off
        
        k = dsearchn(X(:,1),C);
        LogIndexes(ChIndx,:) = sort(X(k,3));
    
    end
end

% To get the Log stamp values. Should actually be placed outside this file
MinIndx = 11;     % min index of log index name to know which index are best for SAGE 
Offset = MinIndx + floor(TimeSize(bestIndx)/2) - 1;
BestLogIndexes = round(mean(LogIndexes,1)) + Offset

PlotIndx = BestLogIndexes - min(MeasChIndxRange);

figure(21)
plot(PhaseCorrectedUnwrapY(1,:))
xline(PlotIndx*200)
ylabel("Phase in radians")
xlabel("Log Sample Index (every sec/skip rest)")
set(gca,"FontSize",14)
title("Unwrapped phase for Channel 0 with stable index positions",'FontSize',18)
    

save('MatFiles/BestLogIndexes.mat','BestLogIndexes');

toc