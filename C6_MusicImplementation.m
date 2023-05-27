%% Contributions by Nikhil Challa
clear;clc; close all;

tic
% Referring to below link to get closest semi-definite matrix
% https://math.stackexchange.com/questions/1098039/converting-a-matrix-to-the-nearest-positive-definite-matrix


StoreAlphaTau = load('MatFiles/StoreAlphaTau3.mat').StoreAlphaTau;
CalibrationResults = load('MatFiles/CalibrationResults.mat').CalibrationResults;

NumTimeSteps = size(StoreAlphaTau,1);
doasFail = zeros(1,NumTimeSteps); % This is to show failure in DOA estimation if we dont do extra steps
doas = zeros(1,NumTimeSteps);
doasSage = zeros(1,NumTimeSteps);
doasCal = zeros(1,NumTimeSteps);
doasSheet = zeros(1,NumTimeSteps);

% We now need to remove the contributions from HW
FinalData = angle(transpose(StoreAlphaTau(:,:,2))) + angle(CalibrationResults);
Indx = FinalData > pi;
FinalData(Indx) = FinalData(Indx) - 2*pi;
Indx = FinalData < -pi;
FinalData(Indx) = FinalData(Indx) + 2*pi;

% Directly from SAGE output
FinalDataSage = angle(transpose(StoreAlphaTau(:,:,2)));

filename = 'Data.xlsx';
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

    CovMatFail = CovMat; % This is to show failure in DOA estimation if we dont do extra steps
    CovMat = F_nearestSPD(CovMat);
    CovMatSage = F_nearestSPD(CovMatSage);
    CovMatCal = F_nearestSPD(CovMatCal);
    CovMatSheet = F_nearestSPD(CovMatSheet);

    [Vfail,Dfail] = eig(CovMat);
    [V,D] = eig(CovMat);
    [Vsage,Dsage] = eig(CovMatSage);
    [Vcal,Dcal] = eig(CovMatCal);
    [Vsheet,Dsheet] = eig(CovMatSheet);


    % We can clearly see the Eigen values has numbers with very small
    % value after decimal point hence rounding helps
    
    % Vfail and Dfail no change
    Dnew = round(D,4);
    DnewSage = round(Dsage,4);
    DnewCal = round(Dcal,4);
    DnewSheet = round(Dsheet,4);


    % We need to introduce some form of resolution reduction to get all DOA
    % points.
    doasFail(timeStep) = F_RecResolReductForMusic(Vfail,Dfail,timeStep);
    doas(timeStep) = F_RecResolReductForMusic(V,Dnew,timeStep);
    doasSage(timeStep) = F_RecResolReductForMusic(Vsage,DnewSage,timeStep);
    doasCal(timeStep) = F_RecResolReductForMusic(Vcal,DnewCal,timeStep);
    doasSheet(timeStep) = F_RecResolReductForMusic(Vsheet,DnewSheet,timeStep);

end

figure(1)
plot(1:NumTimeSteps,Angle,'-.',1:NumTimeSteps,doas,'-x',1:NumTimeSteps,doasFail,'-^',1:NumTimeSteps,doasSage,'-+',1:NumTimeSteps,doasCal,'-o',1:NumTimeSteps,doasSheet,'-*')

xlabel("Position Index")
ylabel("Direction of arrival in Degrees")
legend("Ground Truth","Direction of Arrival Estimation","Direction of Arrival Estimation without Processing","Direction of Arrival Estimation from Sage","Direction of Arrival using Cal data","Direction of Arrival using Sheet data",'Location','best')
yline([0],'-',{"DOA Estimation Failure Line"},'HandleVisibility','off')
set(gca,"FontSize",14)
title("Direction of Arrival plots",'FontSize',18)

toc
