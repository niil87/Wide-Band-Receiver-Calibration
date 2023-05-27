%% Contributions by Nikhil Challa
clear;clc; close all;

tic

%%%%%%%%%%%%%%  INPUTS %%%%%%%%%%%%%%%%%%
K = 1024;    % number of samples
L = 10;      % number of multipath
CIR_RATE = 1; % 
B2BrefIndx = 15;
ChIndx = 0:3;
FirstIndxY = 11;
FirstIndxX = 11;
MeasChIndxRange = 11:635;
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

    disp(["Processing Ch", ChNo])
    X = fft(x(ChNo+1,StartIndx:EndIndx)'); X(1) = 0;

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
                    ylabel("Scaled Amplitude")
                    xlabel("Time Delay")
                    ylim([0 1])
                    set(gca,"FontSize",14)
                    title("Plot of multipath with sequential elimination",'FontSize',18)
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
save('Matfiles/StoreAlphaTau3.mat','StoreAlphaTau');

% incase variables are overwritten.

%% Loading values for plots
StoreAlphaTau = load('Matfiles/StoreAlphaTau3.mat').StoreAlphaTau;


% for unwrapping we use different property as compared to UnwrapPhase.m
% Here we know that phase can only be monotonically increasing before
% midway and monotonically decreasing after midway. Below is mostly visual
% analysis

MidIndx = length(LTSi) / 2;
Indx = 1:20;

figure(31)

for subpIndx = 1:4
    subplot(2,2,subpIndx)
    hold on
    PhaseOrig = angle(StoreAlphaTau(:,subpIndx,2));
    [~,PhaseUnwrap,~] = F_UnwrapPhaseMonotone(PhaseOrig,LTSi,MidIndx);
    plot(Indx,PhaseOrig,'--x',Indx,PhaseUnwrap,'-o')
    xlim([0 21])
    yline([pi,-pi,-3*pi,-5*pi,-7*pi,-9*pi,-11*pi,-13*pi])
    ylabel("Phase in radians")
    xlabel("Stable index")
    legend("Original phase","Unwrapped phase","Location",'best')
    set(gca,"FontSize",14)
    titleStr = "Channel " + string((subpIndx-1));
    title(titleStr,'FontSize',18)
    sgtitle('Phase value after SAGE for LoS Component','FontSize',18) 
    hold off
end


%% Verifying Measurement data (compare it against SAGE output)

filename = 'Data.xlsx';
A = xlsread(filename);

MidIndx = length(LTSi) / 2;
Indx = 1:20;

FinalMeasPhase = zeros(4,length(LTSi));
figure(32)
sgtitle('Phase value from Data sheet','FontSize',18) 
for subpIndx = 1:4
    subplot(2,2,subpIndx)
    hold on
    PhaseOrig = A(2:21,(9 + 4*subpIndx) );
    [~,PhaseUnwrap,~] = F_UnwrapPhaseMonotone(PhaseOrig,LTSi,MidIndx);
    plot(Indx,PhaseOrig,'--x',Indx,PhaseUnwrap,'-o')
    yline([pi,-pi,-3*pi,-5*pi,-7*pi,-9*pi,-11*pi,-13*pi])
    ylabel("Phase in radians")
    xlabel("Stable index")
    legend("Original angle","Angle with unwrap","Location",'best')
    set(gca,"FontSize",14)
    titleStr = "Channel " + string((subpIndx-1));
    title(titleStr,'FontSize',18)
    hold off
    
    FinalMeasPhase(subpIndx,:) = PhaseOrig;

end

save('Matfiles/FinalMeasPhase.mat','FinalMeasPhase');

toc

%% optimizing function

function ArrayMin = MinFn(opt_Tau)

    global K X Y_i
    ArrayMin = -1* abs( ( exp(-1i*2*pi*(0:K-1)'*opt_Tau/K ) .* X )' * Y_i );

end







