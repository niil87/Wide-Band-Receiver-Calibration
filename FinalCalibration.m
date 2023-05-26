%% Contributions by Nikhil Challa
clear;clc; close all;

StoreAlphaTau = load('StoreAlphaTau3.mat').StoreAlphaTau;
LTSi = load('MatFiles/BestLogIndexes.mat').BestLogIndexes;
freq = 5.725*10^9; % carrier frequency

SamplingTime = 2*10^-9;
%SamplingRate = 2*10^-11;
%SamplingRate = 1/200/1024;

ChNoRange = 0:3;  

filename = 'Data_v3.xlsx';
A = xlsread(filename);

% Note we are multiplying the (-ve) phase contribution from propagation to remove
% this phase on the OTA signal, to get phase contribution from the RF
% device (hence the missing sign below)
PhaseCalcCh{1} = exp(A(2:21,13).*1i);
PhaseCalcCh{2} = exp(A(2:21,17).*1i);
PhaseCalcCh{3} = exp(A(2:21,21).*1i);
PhaseCalcCh{4} = exp(A(2:21,25).*1i);
NumPos = length(LTSi);

g = zeros(length(ChNoRange),NumPos);
PhaseUnwrap = zeros(length(ChNoRange),NumPos);

for ChNo = ChNoRange 
    for indx = 1:NumPos
        a_est = StoreAlphaTau(indx,ChNo+1,2);
        a_real = PhaseCalcCh{ChNo+1}(indx);
        g(ChNo+1,indx) = conj(a_est)*a_real/(a_est*conj(a_est));
    end
    PhaseOrig = angle(g(ChNo+1,:));
    %[~,PhaseUnwrap(ChNo+1,:),~] = UnwrapPhaseMonotone(PhaseOrig,LTSi,10);
    [~,PhaseUnwrap(ChNo+1,:),~] = UnwrapPhase(PhaseOrig,1,2);
    
end


figure(1)
subplot(2,2,1)
hold on
plot(LTSi,angle(g(1,:)),'o', LTSi,PhaseUnwrap(1,:),'x')
title("Calibration Info Ch0")
yline([-2*pi,-pi,pi,2*pi],":",{'-2pi','-pi','pi','2pi'},'LabelHorizontalAlignment','left')
legend("Original Angle","Unwrap Angle",'','','','','Location','south')
hold off

subplot(2,2,2)
hold on
plot(LTSi,angle(g(2,:)),'o', LTSi,PhaseUnwrap(2,:),'x')
title("Calibration Info Ch1")
yline([-2*pi,-pi,pi,2*pi],":",{'-2pi','-pi','pi','2pi'},'LabelHorizontalAlignment','left')
legend("Original Angle","Unwrap Angle",'','','','','Location','south')
hold off

subplot(2,2,3)
hold on
plot(LTSi,angle(g(3,:)),'o', LTSi,PhaseUnwrap(3,:),'x')
title("Calibration Info Ch2")
yline([-2*pi,-pi,pi,2*pi],":",{'-2pi','-pi','pi','2pi'},'LabelHorizontalAlignment','left')
legend("Original Angle","Unwrap Angle",'','','','','Location','south')
hold off

subplot(2,2,4)
hold on
plot(LTSi,angle(g(4,:)),'o', LTSi,PhaseUnwrap(4,:),'x')
title("Calibration Info Ch3")
yline([-2*pi,-pi,pi,2*pi],":",{'-2pi','-pi','pi','2pi'},'LabelHorizontalAlignment','left')
legend("Original Angle","Unwrap Angle",'','','','','Location','south')
hold off