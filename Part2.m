clear all
close all
clc

Hs = [7.61 7.51 6.76 6.68 6.25 6.05 5.88 5.70 5.70 5.56];

Tp = [12.5 12.18 11.46 11.49 10.94 10.78 11.1 10.95 10.71 10.28];

g=9.81;

%% 6.

p1 = polyfit(sqrt(Hs/g),Tp,1);

figure()
plot(Hs,Tp)
hold on
plot(Hs, p1(1).*sqrt(Hs./g)+p1(2))
hold off
legend('raw data','polyfit')
xlabel('H_s [m]')
ylabel('T_p [s]')
%legend({'$T_p = f(H_s)$','$T_p = 14.0175 \sqrt{\frac{H_s}{g}} - 0.0396$'}, ... 
    %'Interpreter','latex','Location','Southeast')
enhance_plot('TIMES',16,1.5)

%% 7.

P = 1 - [1:1:10]./20;
% wd=fitdist(Hs','weibul'); %with matlab function
% A = wd.A;
% B = wd.B;

% Hprob=cdf(wd,[5:0.1:8]);


%by hand
X = log(Hs);
Y = log(-log(1-P));

p2 = polyfit(X,Y,1);
b = p2(1);
a = exp(-p2(2)/p2(1));

% figure()
% plot([5:0.1:8],Hprob)
% enhance_plot('TIMES',16,1.5)

figure()
probplot('weibul',Hs)
enhance_plot('TIMES',16,1.5)

%for 50 years
P_50 = 1 - 1/50;
Hs_50 = exp( log(-log(1-P_50))/b + log(a));

figure()
plot(Hs, 1 - exp(-(Hs./a).^b))
xlabel('H_s')
ylabel('P(H_s)')
title('Weibull distribution') 
enhance_plot('TIMES',16,1.5)

%% 8.

Tp_50 = p1(1).*sqrt(Hs_50./g)+p1(2);


