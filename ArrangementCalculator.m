%% Contributions by Nikhil Challa
clear;clc; close all;

sympref('FloatingPointOutput',true)

% all units in meters
T2Rtable = 1.629;   % Minimum distance from Rx table to Tx table
% T2Rtable = 1.8;     % Recommended value to get mathematically friendly numbers

Fc = 5.725 * 10^9;  % Carrier Frequency

% Calculate wavelength
lambda = 3*10^8 / Fc 

% in percentage of lambda
lowerLimit = 0.1;
higherLimit = 0.9;


%% Calculations for lower limit

syms a b

eqn = sqrt(T2Rtable^2 + (a)^2) - T2Rtable == lambda*lowerLimit;

result = solve(eqn);

% you will get an answer of ~ 0.13 meters or 13cm
SlideStep = result(result > 0) 

%% Calculations for upper limit

% SlideStep = 0.1;      % force setting it to 0.1 to check results
% SlideStep = 0.16;     % recommended value for mathematically friendly numbers                     

eqn =  sqrt(b^2 + T2Rtable^2) - sqrt( (b - SlideStep)^2 + T2Rtable^2 ) == lambda*higherLimit;

result = solve(eqn);

tableHalfL = result(result > 0)



%% Calculations for minimum height

FresnelZone = 8.656*sqrt(T2Rtable*10^6/Fc)