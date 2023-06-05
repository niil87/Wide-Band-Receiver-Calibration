%% Contributions by Nikhil Challa
clear;clc; close all;

sympref('FloatingPointOutput',true)

% all units in meters
% T2Rtable = 0.24;   % Minimum distance from Rx table to Tx table
T2Rtable = 2.0;     % Recommended value to get mathematically friendly numbers

Fc = 5.725 * 10^9;  % Carrier Frequency
Nrx = 4;            % Number of Rx receivers
RxS = 2.6/100;       % Separation between adjacent Rx antenna patches

% Calculate wavelength
lambda = 3*10^8 / Fc;

% in percentage of lambda
lowerLimit = 0.1;
higherLimit = 0.9;

%% Rayleigh distance calculator

La = RxS*(Nrx-1)                % Largest separation between Rx antennas
dR = 2*La*La/lambda             % Rayleigh distance


%% Calculations for lower limit

syms a b

eqn = sqrt(T2Rtable^2 + (a)^2) - T2Rtable == lambda*lowerLimit;

result = solve(eqn);

% you will get an answer of ~ 0.15 meters or 15cm
SlideStep = result(result > 0) 

%% Calculations for upper limit

% SlideStep = 0.1;      % force setting it to 0.1 to check results
SlideStep = 0.16;     % recommended value for mathematically friendly numbers                     

eqn =  sqrt(b^2 + T2Rtable^2) - sqrt( (b - SlideStep)^2 + T2Rtable^2 ) == lambda*higherLimit;

result = solve(eqn);

tableHalfL = result(result > 0)



%% Calculations for minimum height

FresnelZone = 8.656*sqrt(T2Rtable*10^6/Fc)