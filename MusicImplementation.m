%% Contributions by Nikhil Challa
clear;clc; close all;

% StoreAlphaTau = load('MatFiles/StoreAlphaTau3.mat').StoreAlphaTau;
% AlphaAfterCorrection = load('MatFiles/AlphaAfterCorrection.mat').AlphaAfterCorrection;

StoreAlphaTau = load('MatFiles/StoreAlphaTau3.mat').StoreAlphaTau;
AlphaAfterCorrection = load('MatFiles/AlphaAfterCorrection.mat').AlphaAfterCorrection;
NumTimeSteps = size(StoreAlphaTau,1);
doas = zeros(1,NumTimeSteps);
doas2 = zeros(1,NumTimeSteps);

filename = 'Data_v3.xlsx';
A = xlsread(filename);
Angle = A(2:21,6);

for timeStep = 1:NumTimeSteps
    Phase = transpose(angle(StoreAlphaTau(timeStep,:,2)));
    Phase2 = angle(AlphaAfterCorrection(:,timeStep));
    A = exp(-1i*Phase);
    A2 = exp(-1i*Phase2);
    CovMat = A*A';
    CovMat2 = A2*A2';
    [V,D] = eig(CovMat)
    [V2,D2] = eig(CovMat2)
    Dnew = D;
    Dnew2 = D2;
    Dnew(Dnew<0.1)=0;
    Dnew2(Dnew2<0.1)=0;
    CovMatNew = V*Dnew*V';
    CovMatNew2 = V2*Dnew2*V2';
    try
        doas(timeStep) = musicdoa(CovMatNew,1);
        disp(["TimeStep:",timeStep,"doas",doas(timeStep)]);
    catch
        disp(["TimeStep:",timeStep," Failed "]);
    end

    try
        doas2(timeStep) = musicdoa(CovMatNew2,1);
        disp(["TimeStep2:",timeStep,"doas2",doas2(timeStep)]);
    catch
        disp(["TimeStep2:",timeStep," Failed "]);
    end

end

figure(1)
plot(1:NumTimeSteps,doas,1:NumTimeSteps,doas2,1:NumTimeSteps,Angle)
legend("Original","WithSlopeAdjustment","Measured Angle",'Location','southwest')