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

[f,a,S] = jonswap(Hs_50, Tp_50, df, fHighCut, gammaJS);

figure()
plot(f,S)
xlabel('Frequency [Hz]')
ylabel('Spectrum density [m^2 s]')
enhance_plot('TIMES',16,1.5)


%% 10.

phi = 2*pi*rand(1,length(f));

[eta,F,M]=FCalc_irregular(f,h,g,D,a,phi);
