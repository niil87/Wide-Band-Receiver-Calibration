%% Contributions by Nikhil Challa
clear;clc; close all;

StoreAlphaTau = load('MatFiles/StoreAlphaTau3.mat').StoreAlphaTau;
AlphaAfterCorrection = load('MatFiles/AlphaAfterCorrection.mat').AlphaAfterCorrection;
NumTimeSteps = size(StoreAlphaTau,1);
doas = zeros(1,NumTimeSteps);


filename = 'Data.xlsx';
A = xlsread(filename);
Angle = A(2:21,7);

for timeStep = 1:NumTimeSteps
    Phase = transpose(angle(StoreAlphaTau(timeStep,:,2)));
    %A = exp(-1i*Phase) + (randn(4,1)/5 + 1i*randn(4,1)/5);
    A = exp(-1i*Phase);
    CovMat = A*A';
    [V,D] = eig(CovMat)
    Dnew = D;
    Dnew(Dnew<0.1)=0;
    CovMatNew = V*Dnew*V'

    try
        doas(timeStep) = musicdoa(CovMatNew,1);
        disp(["TimeStep:",timeStep,"doas",doas(timeStep)]);
    catch
        disp(["TimeStep:",timeStep," Failed "]);
    end


end

figure(1)
plot(1:NumTimeSteps,doas,1:NumTimeSteps,Angle)
legend("Original","WithSlopeAdjustment","Measured Angle",'Location','north')