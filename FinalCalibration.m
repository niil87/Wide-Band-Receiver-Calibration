%% Contributions by Nikhil Challa
clear;clc; close all;

HACK = 0;  % 1 : this hack is to use ch2 for ch0 and ch3 for ch1, 2 : ch2 for all channels

if HACK == 0
    StoreAlphaTau = load('MatFiles/StoreAlphaTau3.mat').StoreAlphaTau;
elseif HACK == 1
    StoreAlphaTau = load('MatFiles/StoreAlphaTau3_HACK1.mat').StoreAlphaTau;
else
    StoreAlphaTau = load('MatFiles/StoreAlphaTau3_HACK2.mat').StoreAlphaTau;
end

FinalMeasPhase = load('Matfiles/FinalMeasPhase.mat','FinalMeasPhase').FinalMeasPhase;
A = xlsread('Data_v5.xlsx');

% 0 : Use slope adjusted phase value from Excel
% 1 : Use Original data or excel sheet directly with no phase adjustment 
USE_ORIG_MEAS_DATA = 0;

freq = 5.725*10^9; % carrier frequency

SamplingTime = 2*10^-9;

ChNoRange = 0:3;  

% Note we are multiplying the (-ve) phase contribution from propagation to remove
% this phase on the OTA signal, to get phase contribution from the RF
% device (hence the missing sign below)
if USE_ORIG_MEAS_DATA == 0
    PhaseCalcCh{1} = exp(FinalMeasPhase(1,:).*1i);
    PhaseCalcCh{2} = exp(FinalMeasPhase(2,:).*1i);
    PhaseCalcCh{3} = exp(FinalMeasPhase(3,:).*1i);
    PhaseCalcCh{4} = exp(FinalMeasPhase(4,:).*1i);
else
    PhaseCalcCh{1} = exp(A(2:21,13).*1i);
    PhaseCalcCh{2} = exp(A(2:21,17).*1i);
    PhaseCalcCh{3} = exp(A(2:21,21).*1i);
    PhaseCalcCh{4} = exp(A(2:21,25).*1i);
end



NumPos = size(StoreAlphaTau,1);

g = zeros(length(ChNoRange),NumPos);
PhaseUnwrap = zeros(length(ChNoRange),NumPos);

for ChNo = ChNoRange 
    for indx = 1:NumPos
        a_est = (StoreAlphaTau(indx,ChNo+1,2));
        a_real = (PhaseCalcCh{ChNo+1}(indx));
        g(ChNo+1,indx) = conj(a_est)*a_real/(a_est*conj(a_est));
    end
    PhaseOrig = angle(g(ChNo+1,:));
    %[~,PhaseUnwrap(ChNo+1,:),~] = UnwrapPhaseMonotone(PhaseOrig,LTSi,10);
    [~,PhaseUnwrap(ChNo+1,:),~] = UnwrapPhase(PhaseOrig,1,2);
    
end

Xindx = 1:NumPos;

figure(1)
for ChNo = ChNoRange 
    subplot(2,2,ChNo+1)
    hold on
    plot(Xindx,angle(g(ChNo+1,:)),'o', Xindx,PhaseUnwrap(ChNo+1,:),'x')
    ylabel("Phase in radians")
    xlabel("Position Index")
    yline([-2*pi,-pi,pi,2*pi],":",{'-2pi','-pi','pi','2pi'},'LabelHorizontalAlignment','left')
    legend("Original Angle","Unwrap Angle",'','','','','Location','best')
    set(gca,"FontSize",14)
    titleStr = "Channel " + string((ChNo)) + " : Calibration Info";
    title(titleStr,'FontSize',18)
    hold off
end


CalibrationResults = g;
save('Matfiles/CalibrationResults.mat','CalibrationResults');

PhaseUnwrap(2,:) = PhaseUnwrap(2,:) + PhaseUnwrap(1,10) - PhaseUnwrap(2,10);
PhaseUnwrap(3,:) = PhaseUnwrap(3,:) + PhaseUnwrap(1,10) - PhaseUnwrap(3,10);
PhaseUnwrap(4,:) = PhaseUnwrap(4,:) + PhaseUnwrap(1,10) - PhaseUnwrap(4,10);

figure(2)
plot(Xindx,PhaseUnwrap(1,:), Xindx,PhaseUnwrap(2,:), Xindx,PhaseUnwrap(3,:), Xindx,PhaseUnwrap(4,:) )
legend("Channel 1","Channel 2","Channel 3","Channel 4",'Location','best')
ylim([-pi pi])
title("Phase Hacked",'FontSize',18)

HackedPhase = exp(1i*PhaseUnwrap);
save('Matfiles/HackedPhase.mat','HackedPhase');