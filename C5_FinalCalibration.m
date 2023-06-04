%% Contributions by Nikhil Challa
clear;clc; close all;

tic

StoreAlphaTau = load('MatFiles/StoreAlphaTau3.mat').StoreAlphaTau;
FinalMeasPhase = load('Matfiles/FinalMeasPhase.mat','FinalMeasPhase').FinalMeasPhase;

freq = 5.725*10^9; % carrier frequency

SamplingTime = 2*10^-9;

ChNoRange = 0:3;  

% Note we are multiplying the (-ve) phase contribution from propagation to remove
% this phase on the OTA signal, to get phase contribution from the RF
% device (hence the missing -ve sign below)
PhaseCalcCh{1} = exp(FinalMeasPhase(1,:).*1i);
PhaseCalcCh{2} = exp(FinalMeasPhase(2,:).*1i);
PhaseCalcCh{3} = exp(FinalMeasPhase(3,:).*1i);
PhaseCalcCh{4} = exp(FinalMeasPhase(4,:).*1i);

NumPos = size(StoreAlphaTau,1);

g = zeros(length(ChNoRange),NumPos);
PhaseUnwrap = zeros(length(ChNoRange),NumPos);

for ChNo = ChNoRange 
    for indx = 1:NumPos
        a_est = (StoreAlphaTau(indx,ChNo+1,2));
        a_meas = (PhaseCalcCh{ChNo+1}(indx));
        g(ChNo+1,indx) = conj(a_est)*a_meas/(a_est*conj(a_est));
    end
    PhaseOrig = angle(g(ChNo+1,:));
    [~,PhaseUnwrap(ChNo+1,:),~] = F_UnwrapPhase(PhaseOrig,1,2);
end

Xindx = 1:NumPos;

figure(51)
tiledlayout(2,2,'Padding','Compact');
sgtitle('Calibration Info','FontSize',26) 
for ChNo = ChNoRange 
    nexttile
    hold on
    plot(Xindx,angle(g(ChNo+1,:)),'o', Xindx,PhaseUnwrap(ChNo+1,:),'x')
    ylabel("Phase in radians")
    xlabel("Position Index")
    yline([-2*pi,-pi,pi,2*pi],":",{'-2pi','-pi','pi','2pi'},'LabelHorizontalAlignment','left')
    legend("Original Angle","Unwrap Angle",'','','','','Location','best')
    set(gca,"FontSize",18)
    titleStr = "Channel " + string((ChNo));
    title(titleStr,'FontSize',22)
    hold off
end


CalibrationResults = g;
save('Matfiles/CalibrationResults.mat','CalibrationResults');

toc