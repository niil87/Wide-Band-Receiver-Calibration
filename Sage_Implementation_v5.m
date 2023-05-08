%% Contributions by Nikhil Challa
clear;clc; close all;

%%%%%%%%%%%%%%  INPUTS %%%%%%%%%%%%%%%%%%
K = 1024;    % number of samples
L = 10;      % number of multipath
CIR_RATE = 1; % 
B2BrefIndx = 15;
ChIndx = 0:3;
FirstIndxY = 11;
FirstIndxX = 11;
MeasChIndxRange = 11:635;
HACK = 0;  % this hack is to use ch2 for ch0 and ch3 for ch1
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global K X Y_i

format long

LTSi = load('MatFiles/BestLogIndexes.mat').BestLogIndexes;

%% Main Code

% (number of meas points, number of channels, number of measurements per channel)
StoreAlphaTau = zeros(length(LTSi),4,2);
ChPhase = zeros(1,length(LTSi));
ChAmp = zeros(1,length(LTSi));


%% Read frequency corrected samples
x = load('MatFiles/x_fcorr.mat').x_fcorr;
y = load('MatFiles/y_fcorr.mat').y_fcorr;


MaxIndx = 0;
for ChNo = ChIndx

    StartIndx = (B2BrefIndx - FirstIndxX)*K*200 + 1;
    EndIndx = StartIndx + K - 1;

    if ((ChNo < 2) && (HACK == 1))
        disp(["Doing hack, using Ch", ChNo+2, " for Ch", ChNo])
        X = fft(x(ChNo+3,StartIndx:EndIndx)'); X(1) = 0;
    else
        disp(["Processing Ch", ChNo])
        X = fft(x(ChNo+1,StartIndx:EndIndx)'); X(1) = 0;
    end

    options = optimoptions('fminunc',"Display","off");

    for SrNo = 1:length(LTSi)

            fprintf('SrNo %d was processed at time %s\n',SrNo, datestr(now,'HH:MM:SS.FFF'))
            StartIndx = (LTSi(SrNo) - FirstIndxY)*K*200 + 1;
            EndIndx = StartIndx + K - 1;
            Y = fft(y(ChNo+1, StartIndx : EndIndx )'); Y(1) = 0;
            
            % compute multipath and Implemeting SAGE
            plotLimit = K/8;
            TauList = (0:0.2:K - 0.2);
            
            % initializing complex amplitude shift to low complex number
            Alpha_hat =  zeros(1,L);
            Tau_hat = zeros(1,L);
            itrLoops = 5;             % number of iterations
            Tau_init = zeros(1,L);
            lockIndex = 0;
            
            m = (0:1:K-1);
            A = X.' .* exp(-1i*2*pi/K*m'*m);
            Ac = conj(A);
            
            for iter = 1:itrLoops
                % Looping over all Multipath for estimation
                for i = 1 : L
            
                    if iter == 1
                        FvalMin = -10^10;
                        Bestindx = 0;
                        for j = 1:length(TauList)
                            XsumOpt = zeros(K,1);
                            for Mpc = 1 : L 
                                if Mpc ~= i
                                    XsumOpt = XsumOpt + Alpha_hat(Mpc) .* ( exp( -1i*2*pi*(0:K-1)'*Tau_hat(Mpc)/K ) .* X );
                                end
                            end
                            Y_i = Y - XsumOpt;
                           
                            Val = abs( ( exp(-1i*2*pi*(0:K-1)'*TauList(j)/K ) .* X )' * Y_i );
                            
                            if Val > FvalMin
                                FvalMin = Val;
                                Tau_hat(i) = TauList(j);
                                Tau_init(i) = TauList(j);
                                lockIndex = j;
                            end
                        end
                    end

                    XsumOpt = zeros(K,1);
                    for Mpc = 1 : L 
                        if Mpc ~= i
                            %disp(["facotring in:" , Mpc]);
                            XsumOpt = XsumOpt + Alpha_hat(Mpc) .* ( exp( -1i*2*pi*(0:K-1)'*Tau_hat(Mpc)/K ) .* X );
                        end
                    end
                    Y_i = Y - XsumOpt;
            
            
                    %%%%%%%%%%%%%% plotting after subracting %%%%%%%%%%%%
                    Pwr = abs( Ac*Y_i/(norm(X,"fro")^2) ) ;
                    figure(1)
                    plot(Pwr(1:plotLimit))
                    title("Plot of multipath with sequential elimination")
                    grid on
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
                    [Tau_hat(i),fval] = fminunc(@MinFn,Tau_hat(i),options);
            
                    X_i = X .* exp( -1i*2*pi*(0:K-1)' * Tau_hat(i) / K );
                    Alpha_hat(i) = X_i' * Y_i / norm(X_i,"fro")^2;
                    
                end
                [TauSort,Indx] = sort(Tau_hat);
            end
            [Val, MaxIndx] = max(abs(Alpha_hat));
            StoreAlphaTau(SrNo,ChNo+1,1) = Tau_hat(MaxIndx);
            StoreAlphaTau(SrNo,ChNo+1,2) = Alpha_hat(MaxIndx);

    end
end

% Storing calculated values
if HACK == 1
    save('Matfiles/StoreAlphaTau3_HACK.mat','StoreAlphaTau');
else
    save('Matfiles/StoreAlphaTau3.mat','StoreAlphaTau');
end

% incase variables are overwritten.

%% Loading values for plots
if HACK == 1
    StoreAlphaTau = load('Matfiles/StoreAlphaTau3_HACK.mat').StoreAlphaTau;
else
    StoreAlphaTau = load('Matfiles/StoreAlphaTau3.mat').StoreAlphaTau;
end

% for unwrapping we use different property as compared to UnwrapPhase.m
% Here we know that phase can only be monotonically increasing before
% midway and monotonically decreasing after midway. Below is mostly visual
% analysis

MidIndx = length(LTSi) / 2;

figure(31)
subplot(2,2,1)
hold on
PhaseOrig = angle(StoreAlphaTau(:,1,2));
[~,PhaseUnwrap,~] = UnwrapPhaseMonotone(PhaseOrig,LTSi,MidIndx);
plot(LTSi,PhaseOrig,LTSi,PhaseUnwrap,'o')
xline(LTSi)
yline([-pi,pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi])
ylabel("Phase in radians")
xlabel("Samples/CIR-RATE")
legend("Original angle","Angle with unwrap")
title("Channel 0 : Phase value after SAGE for LoS Component")
hold off

subplot(2,2,2)
hold on
PhaseOrig = angle(StoreAlphaTau(:,2,2));
[~,PhaseUnwrap,~] = UnwrapPhaseMonotone(PhaseOrig,LTSi,MidIndx);
plot(LTSi,PhaseOrig,LTSi,PhaseUnwrap,'o')
xline(LTSi)
yline([-pi,pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi])
ylabel("Phase in radians")
xlabel("Samples/CIR-RATE")
legend("Original angle","Angle with unwrap")
title("Channel 1 : Phase value after SAGE for LoS Component")
hold off

subplot(2,2,3)
hold on
PhaseOrig = angle(StoreAlphaTau(:,3,2));
[~,PhaseUnwrap,~] = UnwrapPhaseMonotone(PhaseOrig,LTSi,MidIndx);
plot(LTSi,PhaseOrig,LTSi,PhaseUnwrap,'o')
xline(LTSi)
yline([-pi,pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi])
ylabel("Phase in radians")
xlabel("Samples/CIR-RATE")
legend("Original angle","Angle with unwrap")
title("Channel 2 : Phase value after SAGE for LoS Component")
hold off

subplot(2,2,4)
hold on
PhaseOrig = angle(StoreAlphaTau(:,4,2));
[~,PhaseUnwrap,~] = UnwrapPhaseMonotone(PhaseOrig,LTSi,MidIndx);
plot(LTSi,PhaseOrig,LTSi,PhaseUnwrap,'o')
xline(LTSi)
yline([-pi,pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi])
ylabel("Phase in radians")
xlabel("Samples/CIR-RATE")
legend("Original angle","Angle with unwrap")
title("Channel 3 : Phase value after SAGE for LoS Component")
hold off

%% Verifying Measurement data against SAGE output

filename = 'Data_v4.xlsx';
A = xlsread(filename);

MidIndx = length(LTSi) / 2;

figure(32)
subplot(2,2,1)
hold on
PhaseOrig = A(2:21,13);
[~,PhaseUnwrap,~] = UnwrapPhaseMonotone(PhaseOrig,LTSi,MidIndx);
plot(LTSi,PhaseOrig,LTSi,PhaseUnwrap,'o')
xline(LTSi)
yline([-pi,pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi])
ylabel("Phase in radians")
xlabel("Samples/CIR-RATE")
legend("Original angle","Angle with unwrap")
title("Channel 0 : Phase value from Data sheet")
hold off

subplot(2,2,2)
hold on
PhaseOrig = A(2:21,17);
[~,PhaseUnwrap,~] = UnwrapPhaseMonotone(PhaseOrig,LTSi,MidIndx);
plot(LTSi,PhaseOrig,LTSi,PhaseUnwrap,'o')
xline(LTSi)
yline([-pi,pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi])
ylabel("Phase in radians")
xlabel("Samples/CIR-RATE")
legend("Original angle","Angle with unwrap")
title("Channel 1 : Phase value from Data sheet")
hold off

subplot(2,2,3)
hold on
PhaseOrig = A(2:21,21);
[~,PhaseUnwrap,~] = UnwrapPhaseMonotone(PhaseOrig,LTSi,MidIndx);
plot(LTSi,PhaseOrig,LTSi,PhaseUnwrap,'o')
xline(LTSi)
yline([-pi,pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi])
ylabel("Phase in radians")
xlabel("Samples/CIR-RATE")
legend("Original angle","Angle with unwrap")
title("Channel 1 : Phase value from Data sheet")
hold off

subplot(2,2,4)
hold on
PhaseOrig = A(2:21,25);
[~,PhaseUnwrap,~] = UnwrapPhaseMonotone(PhaseOrig,LTSi,MidIndx);
plot(LTSi,PhaseOrig,LTSi,PhaseUnwrap,'o')
xline(LTSi)
yline([-pi,pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi])
ylabel("Phase in radians")
xlabel("Samples/CIR-RATE")
legend("Original angle","Angle with unwrap")
title("Channel 1 : Phase value from Data sheet")
hold off


%% optimizing function

function ArrayMin = MinFn(opt_Tau)

    global K X Y_i
    ArrayMin = -1* abs( ( exp(-1i*2*pi*(0:K-1)'*opt_Tau/K ) .* X )' * Y_i );

end







