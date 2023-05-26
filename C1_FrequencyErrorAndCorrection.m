%% Contributions by Nikhil Challa
clear;clc; close all;

tic

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
B2BrefIndxRange = 11:20;    %% specify the B2B file name index range eg. cir_ch0_xx.dat
ChNoRange = 0:3;            %% specify all the channels you are intersted in. 0-3 are valid values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumChannels = length(ChNoRange);
SizeB2Bref = (max(B2BrefIndxRange) - min(B2BrefIndxRange) + 1);
K = 1024;            % number of subcarriers per time stamp
CIR_RATE = 200;
format long

mkdir MatFiles;

%% identifying frequency offset and applying frequency correction

PhaseOffsetX = zeros(NumChannels,CIR_RATE*SizeB2Bref);
PhaseOffsetOrigX = zeros(NumChannels,CIR_RATE*SizeB2Bref);
PhaseCorrectedX = zeros(NumChannels,CIR_RATE*SizeB2Bref);
x_fcorr = zeros(NumChannels,CIR_RATE*K*SizeB2Bref);
freqOffset = zeros(1,4);

for ChNo = ChNoRange
    x = zeros(CIR_RATE*K*SizeB2Bref,1);
    ConstPhaseShift = 0;
    for B2BrefIndx = B2BrefIndxRange

        B2BStr = strcat("meas_nikhil_calb_ch",string(ChNo),"/meas_nikhil_calb_ch",string(ChNo),"/cir_ch",string(ChNo), "_", string(B2BrefIndx), ".dat");
        fid2 = fopen(B2BStr);

        raw = fread(fid2, CIR_RATE*K*2,'int16');
        fclose(fid2);

        I = raw(1:2:end)/(2^15-1);
        Q = raw(2:2:end)/(2^15-1);
        
        StrIndx = (B2BrefIndx - min(B2BrefIndxRange))*CIR_RATE*K + 1;
        EndIndx = StrIndx + CIR_RATE*K - 1;

        x( StrIndx : EndIndx ) = I + 1i * Q;

    end

    [PhaseOffsetOrigX(ChNo+1,:),PhaseOffsetX(ChNo+1,:),~] = F_UnwrapPhase(x,K,1);

    % figure(1 + 2*ChNo)
    figure(1)
    subplot(2,2,1 + ChNo)
    plot(PhaseOffsetX(ChNo+1,:))
    ylabel("Phase in radians")
    xlabel("Log Sample Index (sample time)")
    ylim([0 50])
    set(gca,"FontSize",14)
    titleStr = "Channel No:" + string(ChNo);
    title(titleStr,'FontSize',18)
    sgtitle('Phase value before frequency error correction - B2B Data','FontSize',18)  
    
    c = polyfit(0:CIR_RATE*SizeB2Bref - 1,PhaseOffsetX(ChNo+1,:),1)
    
    freqOffset(ChNo+1) = c(1)/(2*pi*1/CIR_RATE);

    % Frequency correction
    for i = 0:CIR_RATE*SizeB2Bref - 1
        StrIndx = i*K + 1;
        EndIndx = (i+1)*K - 1;
        x_fcorr(ChNo + 1, StrIndx:EndIndx,1) = x(StrIndx:EndIndx)*exp(-1i*2*pi*freqOffset(ChNo+1)*i/CIR_RATE);
    end

    [~,PhaseCorrectedX(ChNo+1,:),~] = F_UnwrapPhase(x_fcorr(ChNo + 1,:),K,1);
    
    %figure(2 + 2*ChNo)
    figure(2)
    subplot(2,2,1 + ChNo)
    plot(PhaseCorrectedX(ChNo+1,:))
    ylabel("Phase in radians")
    xlabel("Log Sample Index (sample time)")
    ylim([-0.02 0.02])
    set(gca,"FontSize",14)
    titleStr = "Channel No:" + string(ChNo);
    title(titleStr,'FontSize',18)
    sgtitle('Phase value after frequency error correction - B2B Data','FontSize',18) 

end

save('MatFiles/freqOffset.mat','freqOffset');
save('MatFiles/PhaseOffsetOrigX.mat','PhaseOffsetOrigX');
save('MatFiles/PhaseCorrectedX.mat','PhaseCorrectedX');
save('MatFiles/x_fcorr.mat','x_fcorr');

%% Performing Frequency error correction For Y 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INPUTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MeasChIndxRange = 11:635;    %% specify the B2B file name index range eg. cir_ch0_xx.dat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SizeB2Bref = (max(MeasChIndxRange) - min(MeasChIndxRange) + 1);
PhaseCorrectedY = zeros(NumChannels,CIR_RATE*SizeB2Bref);
PhaseCorrectedUnwrapY = zeros(NumChannels,CIR_RATE*SizeB2Bref);
CorrPhase = zeros(NumChannels,CIR_RATE*SizeB2Bref);
y_fcorr = zeros(NumChannels,SizeB2Bref*CIR_RATE*K);

for ChNo = ChNoRange
    
    ConstPhaseShift = 0;

    disp(["Working on channel:",ChNo]);
    for refIndx = MeasChIndxRange

        Str = strcat("meas_nikhil_array/meas_nikhil_array/cir_ch",string(ChNo), "_", string(refIndx), ".dat");
        fid2 = fopen(Str);

        raw = fread(fid2, CIR_RATE*K*2,'int16');
        fclose(fid2);
        
        I = raw(1:2:end)/(2^15-1);
        Q = raw(2:2:end)/(2^15-1);

        y_temp = I + 1i * Q;

        y = zeros(CIR_RATE*K,1);

        % Frequency correction
        for i = 0:(CIR_RATE-1)
            StrIndx = i*K + 1;
            EndIndx = StrIndx + K - 1;
            y(StrIndx:EndIndx,1) = y_temp(StrIndx:EndIndx)*exp(-1i*2*pi*freqOffset(ChNo+1)*((refIndx - min(MeasChIndxRange))*CIR_RATE + i)/CIR_RATE);
        end

        StrIndx = (refIndx - min(MeasChIndxRange))*CIR_RATE*K + 1;
        EndIndx = StrIndx + CIR_RATE*K - 1;
        y_fcorr(ChNo + 1,StrIndx : EndIndx ) = y;

    end 

    [PhaseCorrectedY(ChNo+1,:),PhaseCorrectedUnwrapY(ChNo+1,:),CorrPhase(ChNo+1,:)] = F_UnwrapPhase(y_fcorr(ChNo + 1,:),K,1);

    figure(3 + ChNo)
    hold on
    plot(MeasChIndxRange,PhaseCorrectedY(ChNo+1,1:CIR_RATE:end),MeasChIndxRange,PhaseCorrectedUnwrapY(ChNo+1,1:CIR_RATE:end))
    yline([-pi,pi,2*pi,3*pi,4*pi,5*pi,6*pi,7*pi,8*pi,9*pi])
    xlabel("Log Sample Index (every sec/skip rest)")
    ylabel("Phase in radians")
    legend("Original Phase","Unwrapped phase")
    set(gca,"FontSize",14)
    titleStr = "Channel No:" + string(ChNo) + "  Phase value after removing frequency offset - OTA Data";
    title(titleStr,'FontSize',18)
    hold off

end

save('MatFiles/y_fcorr.mat','y_fcorr','-v7.3');
save('MatFiles/PhaseCorrectedY.mat','PhaseCorrectedY');
save('MatFiles/PhaseCorrectedUnwrapY.mat','PhaseCorrectedUnwrapY');
save('MatFiles/CorrPhase.mat','CorrPhase');

toc
