clear all
close all
clc

h = 30; %depth m
D = 6; %pile diameter m
g = 9.81; %gravity m/s^2

%% Part 1
%% 2.
H = 6; %height m
T = 12; %period s

omega = 2*pi/T;
f = 1/T;
k = kSolve(f,g,h);
x=0;
L = 2*pi/k;

t = linspace(0,100,500);

eta = H/2*cos(omega.*t - k*x);
u_z0 = omega.*H/2*cosh(k*h)/sinh(k*h)*cos(omega.*t); %horizontal velocity
dudt_z0 = -omega.^2.*H/2*cosh(k*h)/sinh(k*h)*sin(omega.*t); %horizontal acc




figure()
hold on
plot(t,eta)
plot(t,u_z0)
plot(t,dudt_z0)
hold off
xlabel('time')
%ylabel('','Fontsize',14)
legend({'$\eta$ (t)','$u(t)_{|z=0}$','$\frac{\partial u}{\partial t}_{|z=0}$'}, ... 
    'Interpreter','latex','Location','Southeast')
title('Time series')
enhance_plot('TIMES',16,1.5)

%eta et u_z0 même phase, dudt_z0 déphasage
%eta augmente, dp augmente, mouvement -> vitesse augmente
%quand eta = 0 acceleration purement verticale et maximale donc u =0 quand eta
% max, acceleration purement horizontale et decroissante donc u=max
% schéma trajectoire particule


%% 3.
%D/L < 0.2 so we can use Morison equation

a = H/2;
t = [1:200];
[F_drag_w,F_inert_w,F_tot_w] = FCalc(f,h,g,D,a,true,40);


[F_drag,F_inert,F_tot] = FCalc(f,h,g,D,a,false,40);

dif = abs(F_tot - F_tot_w)/F_tot; %e^-4
%explain
% H<<h error due to exp not big enough

figure()
plot(t,10^(-6)*F_drag,t,10^(-6)*F_inert)
xlabel('time','Fontsize',14)
ylabel('Inline force [MN]','Fontsize',14)
legend({'Drag force','Inertial force'},'Fontsize',12)
title('Drag and inertial forces','Fontsize',14)
enhance_plot('TIMES',16,1.5)

figure()
plot(t,10^(-6)*F_tot,'g',t,10^(-6)*F_tot_w,'k--')
xlabel('time','Fontsize',14)
ylabel('Inline force [MN]','Fontsize',14)
legend({'Without wheeler stretching','With wheeler stretching'},'Fontsize',12)
title('Total inline force','Fontsize',14)
enhance_plot('TIMES',16,1.5)

figure()
plot(t,10^(-3)*F_drag,'g',t,10^(-3)*F_drag_w,'k--')
xlabel('time','Fontsize',14)
ylabel('Force [kN]','Fontsize',14)
legend({'Without wheeler stretching','With wheeler stretching'},'Fontsize',12)
title('Drag force','Fontsize',14)
enhance_plot('TIMES',16,1.5)

figure()
plot(t,10^(-6)*F_inert,'g',t,10^(-6)*F_inert_w,'k--')
xlabel('time','Fontsize',14)
ylabel('Force [MN]','Fontsize',14)
legend({'Without wheeler stretching','With wheeler stretching'},'Fontsize',12)
title('Inertial force','Fontsize',14)
enhance_plot('TIMES',16,1.5)

%% 4.
%according to tempel et machin with D/L=0.04 CM ~2 

%% 5 
eta = a*cos(omega.*t);

[~,~,F_tot,M] = FCalc(f,h,g,D,a,false,40);

figure()
plot(t,eta,t,10^(-6)*F_tot,t,10^(-6)*M/2)
xlabel('time','Fontsize',14)
ylabel('','Fontsize',14)
legend({'Free surface elevation','Inline force','Overturning moment'},'Fontsize',12)
title('Time series','Fontsize',14)
enhance_plot('TIMES',16,1.5)

%F and M 

NUM = [10 15 20 25 30 35 40 50 70 90]; %trying different number of points in the horizontal direction


for i =1:(length(NUM))
    [~,~,F_tot5(i,:),M5(i,:)] = FCalc(f,h,g,D,a,false,NUM(i));
end

figure()

subplot(1,2,1)
hold on
plot(NUM, 10^(-6)*max(F_tot5'))
xlabel('Points along the vertical')
ylabel('F [MN]')
title('Inline force')
enhance_plot('TIMES',16,1.5)
plot(NUM,10^(-6)*max(F_tot5(end,:))*ones(1,length(NUM)),'k')
hold off

subplot(1,2,2)
plot(NUM, 10^(-6)*max(M5'))
xlabel('Points along the vertical')
ylabel('M [MN m]') %check unit
title('Overturning moment')
enhance_plot('TIMES',16,1.5)
hold on
plot(NUM,10^(-6)*max(M5(end,:))*ones(1,length(NUM)),'k')


