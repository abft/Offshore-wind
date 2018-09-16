clear all
close all
clc

h = 30;
g = 9.81;
D = 6;

%% 9.
load('data_part2.mat')

T_dur = 3600; %1h time serie
df = 1/T_dur;
fHighCut = 0.5; %1/2dt with dt =1
gammaJS = 3.3; %peak enhancement factor
t = [1:1:T_dur];

[f,a,S] = jonswap(Hs_50, Tp_50, df, fHighCut, gammaJS);

figure()
plot(f,S)
xlabel('Frequency [Hz]')
ylabel('Spectrum density [m^2 s]')
title('Spectrum function')
enhance_plot('TIMES',16,1.5)


%% 10.

phi = 2*pi*rand(1,length(f));

%[eta,F,M]=FCalc_irregular(f,h,g,D,a,phi);

%save('data_part3','eta','F','M');
load('data_part3')

figure()
plot(t,eta)
xlabel('time [s]')
ylabel('Surface elevation [m]')
xlim([0,3600])
title('\eta')
enhance_plot('TIMES',16,1.5)

figure()
plot(t,F)
xlabel('time [s]')
ylabel('Inline force [N]')
xlim([0,3600])
title('F')
enhance_plot('TIMES',16,1.5)

figure()
plot(t,M)
xlabel('time [s]')
ylabel('Overturning moment [Nm]')
xlim([0,3600])
title('M')
enhance_plot('TIMES',16,1.5)

STD = [std(eta),std(F),std(M)];
MEAN = [mean(eta),mean(F),mean(M)];
MAX= [max(eta),max(F),max(M)];
MIN = [min(eta),min(F),min(M)]; 

RATIO = MAX./STD;
