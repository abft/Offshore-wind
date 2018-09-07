clear all
close all
clc

h = 30; %depth m
D = 6; %pile diameter m
g = 9.81; %gravity m/s^2

%% Part 1
% 2.
H = 6; %height m
T = 12; %period s

omega = 2*pi/T;
f = 1/T;
k = kSolve(f,g,h);
x=0;

t = linspace(0,100,1000);

eta = H/2*cos(omega.*t - k*x);
u_z0 = omega.*H/2*cosh(k*h)/sinh(k*h)*cos(omega.*t); %horizontal velocity
dudt_z0 = -omega.^2.*H/2*cosh(k*h)/sinh(k*h)*sin(omega.*t); %horizontal acc

figure()
hold on
plot(t,eta)
plot(t,u_z0)
plot(t,dudt_z0)
hold off
xlabel('time','Fontsize',14)
%ylabel('','Fontsize',14)
legend({'$\eta$ (t)','$u(t)_{|z=0}$','$\frac{\partial u}{\partial t}_{|z=0}$'}, ... 
    'Interpreter','latex','Location','Southeast','Fontsize',14)
title('Time series','Fontsize',14)

%eta et u_z0 même phase, dudt_z0 déphasage
%eta augmente, dp augmente, mouvement -> vitesse augmente
%quand eta = 0 acceleration purement verticale et maximale donc u =0 quand eta
% max, acceleration purement horizontale et decroissante donc u=max
% schéma trajectoire particule



