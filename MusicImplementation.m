%% Contributions by Nikhil Challa
clear;clc; close all;

% Referring to below link to get closest semi-definite matrix
% https://math.stackexchange.com/questions/1098039/converting-a-matrix-to-the-nearest-positive-definite-matrix

HARD_FAIL_ENABLE = 0;
USE_SAGE_DIRECT = 0;
USE_HACKED_CAL_DATA = 0;

HACK = 0;  % 1 : this hack is to use ch2 for ch0 and ch3 for ch1, 2 : ch2 for all channels

% Importing values
if HACK == 0
    StoreAlphaTau = load('MatFiles/StoreAlphaTau3.mat').StoreAlphaTau;
elseif HACK == 1
    StoreAlphaTau = load('MatFiles/StoreAlphaTau3_HACK1.mat').StoreAlphaTau;
else
    StoreAlphaTau = load('MatFiles/StoreAlphaTau3_HACK2.mat').StoreAlphaTau;
end

if USE_HACKED_CAL_DATA == 0
    CalibrationResults = load('MatFiles/CalibrationResults.mat').CalibrationResults;
else
    CalibrationResults = load('Matfiles/HackedPhase.mat','HackedPhase').HackedPhase;
end

NumTimeSteps = size(StoreAlphaTau,1);
doas = zeros(1,NumTimeSteps);
doasSage = zeros(1,NumTimeSteps);
doasCal = zeros(1,NumTimeSteps);
doasSheet = zeros(1,NumTimeSteps);

% We now need to remove the contributions from HW
FinalData = angle(transpose(StoreAlphaTau(:,:,2))) - angle(CalibrationResults);
Indx = FinalData > pi;
FinalData(Indx) = FinalData(Indx) - 2*pi;
Indx = FinalData < -pi;
FinalData(Indx) = FinalData(Indx) + 2*pi;

% Directly from SAGE output
FinalDataSage = angle(transpose(StoreAlphaTau(:,:,2)));


filename = 'Data_v5.xlsx';
A = xlsread(filename);
Angle = A(2:21,7);

PhaseFromSheet(1,:) = A( 2:21, (9 + 4*1) );
PhaseFromSheet(2,:) = A( 2:21, (9 + 4*2) );
PhaseFromSheet(3,:) = A( 2:21, (9 + 4*3) );
PhaseFromSheet(4,:) = A( 2:21, (9 + 4*4) );

for timeStep = 1:NumTimeSteps
    Phase = FinalData(:,timeStep);
    PhaseSage = FinalDataSage(:,timeStep);
    PhaseCal = angle(CalibrationResults(:,timeStep));
    PhaseSheet = PhaseFromSheet(:,timeStep);

    B = exp(-1i*Phase);
    Bsage = exp(-1i*PhaseSage);
    BCal = exp(-1i*PhaseCal);
    Bsheet = exp(-1i*PhaseCal);

    CovMat = B*B';
    CovMatSage = Bsage*Bsage';
    CovMatCal = BCal*BCal';
    CovMatSheet = Bsheet*Bsheet';

    CovMat = nearestSPD(CovMat);
    CovMatSage = nearestSPD(CovMatSage);
    CovMatCal = nearestSPD(CovMatCal);
    CovMatSheet = nearestSPD(CovMatSheet);

    [V,D] = eig(CovMat);
    [Vsage,Dsage] = eig(CovMatSage);
    [Vcal,Dcal] = eig(CovMatCal);
    [Vsheet,Dsheet] = eig(CovMatSheet);

    Dnew = D;
    DnewSage = Dsage;
    DnewCal = Dcal;
    DnewSheet = Dsheet;

    Dnew = round(Dnew,4);
    DnewSage = round(DnewSage,4);
    DnewCal = round(DnewCal,4);
    DnewSheet = round(DnewSheet,4);

    CovMatNew = V*Dnew*V';
    CovMatNewSage = Vsage*DnewSage*Vsage';
    CovMatNewCal = Vcal*DnewCal*Vcal';
    CovMatNewSheet = Vsheet*DnewSheet*Vsheet';

    if HARD_FAIL_ENABLE == 0
        try
            doas(timeStep) = musicdoa(CovMatNew,1);
            disp(["TimeStep:",timeStep,"doasSage",doas(timeStep)]);
        catch
            disp(["TimeStep:",timeStep," Failed "]);
        end

        try
            doasSage(timeStep) = musicdoa(CovMatNewSage,1);
            disp(["TimeStep:",timeStep,"doasSage",doasSage(timeStep)]);
        catch
            disp(["TimeStep:",timeStep," Failed "]);
        end

        try
            doasCal(timeStep) = musicdoa(CovMatNewCal,1);
            disp(["TimeStep:",timeStep,"doasSage",doasCal(timeStep)]);
        catch
            disp(["TimeStep:",timeStep," Failed "]);
        end

        try
            doasSheet(timeStep) = musicdoa(CovMatNewSheet,1);
            disp(["TimeStep:",timeStep,"doasSage",doasSheet(timeStep)]);
        catch
            disp(["TimeStep:",timeStep," Failed "]);
        end

    else
        doas(timeStep) = musicdoa(CovMatNew,1);
        doasSage(timeStep) = musicdoa(CovMatNewSage,1);
        doasCal(timeStep) = musicdoa(CovMatNewCal,1);
        doasSheet(timeStep) = musicdoa(CovMatNewSheet,1);

        disp(["TimeStep:",timeStep,"doas",doas(timeStep)]);
        disp(["TimeStep:",timeStep,"doasSage",doasSage(timeStep)]);
        disp(["TimeStep:",timeStep,"doasSage",doasCal(timeStep)]);
        disp(["TimeStep:",timeStep,"doasSage",doasSheet(timeStep)]);
    end

end


figure(1)
plot(1:NumTimeSteps,doas,1:NumTimeSteps,doasSage,1:NumTimeSteps,doasCal,1:NumTimeSteps,doasSheet,1:NumTimeSteps,Angle)

xlabel("Position Index")
ylabel("Direction of arrival in Degrees")
legend("Direction of Arrival Estimation","Direction of Arrival Estimation from Sage","Direction of Arrival using Cal data","Direction of Arrival using Sheet data","Ground Truth",'Location','best')
set(gca,"FontSize",14)
title("Direction of Arrival plots",'FontSize',18)

figure(2)

subplot(4,1,1)
hold on
plot(angle(transpose(StoreAlphaTau(:,1,2))))
plot(angle(CalibrationResults(1,:)))
plot(FinalData(1,:))
ylim([-pi pi])
hold off

subplot(4,1,2)
hold on
plot(angle(transpose(StoreAlphaTau(:,2,2))))
plot(angle(CalibrationResults(2,:)))
plot(FinalData(2,:))
ylim([-pi pi])
hold off

subplot(4,1,3)
hold on
plot(angle(transpose(StoreAlphaTau(:,3,2))))
plot(angle(CalibrationResults(3,:)))
plot(FinalData(3,:))
ylim([-pi pi])
hold off

subplot(4,1,4)
hold on
plot(angle(transpose(StoreAlphaTau(:,4,2))))
plot(angle(CalibrationResults(4,:)))
plot(FinalData(4,:))
ylim([-pi pi])
hold off