% Report 01 : WaveLoads Monopile
% Part 1 : Linear monochromatic waves

%% Initiliasation of the workspace
clear all
close all
clc

%% Initilisation of the global parameters.
h = 30; % depth (m)
D = 6; % monopile diameter (m)
g = 9.81; % gravity (m/s^2)
H = 6; % Regular wave height (m)
T = 12; % Regular wave period (s)
f = 1/T; % wave frequency (Hz)

%% Question 1

k = kSolve(f,g,h); % calcul of the wavenumber thanks to kSolve.m
L = 2*pi/k; % Wavelength

%% Question 2: Plot time series of eta, u at z=0 and dudt at z=0.

omega = 2*pi/T; % wave angular frequency 
f = 1/T; % wave frequency

x=0; % We look the time series at x=0 (would just be a different phase compare to other x)
t = linspace(0,100,500); % initialisation of the time vector

% Calculate eta, u_z0, dudt_z0 thanks to the linear theory (Aire waves)
eta = H/2*cos(omega.*t-k*x); % surface elevation
u_z0 = omega.*H/2*cosh(k*h)/sinh(k*h)*cos(omega.*t-k*x); % horizontal velocity
dudt_z0 = -omega.^2.*H/2*cosh(k*h)/sinh(k*h)*sin(omega.*t-k*x); %horizontal acc

% Compute Keulegan Carpenter nulber

Um=max(u_z0);
Kc=Um*T/D;

% Plot the time series
fQ2=figure();
plot(t,eta,t,u_z0,t,dudt_z0)
xlabel('t(s)')
legend({'$\eta$ [m]','$u(t)_{|z=0}$ [m/s]','$\frac{\partial u}{\partial t}_{|z=0} [m/s^2]$ '}, ... 
    'Interpreter','latex','Location','Southeast')
title('Time series')
enhance_plot('TIMES',16,1.5)
saveas(fQ2,'Q2.png');



%% Question 3.
% D/L < 0.2 so we can use Morison equation (slender body)

a = H/2; % wave amplitude
t_Fcalc=[1:1:200]; % time vector used in FCalc

% Calcul the inline forces with Wheeler stretching
[F_drag_w,F_inert_w,F_tot_w,M_w] = FCalc(f,h,g,D,a,true,40);

% Calcul the inline forces without Wheeler stretching
[F_drag,F_inert,F_tot,M] = FCalc(f,h,g,D,a,false,40);

% Quantify the Wheeler stretching effect 
Wheeler_effect = abs(max(F_tot) - max(F_tot_w))/max(F_tot); 

% Plot

% Drag and inertial forces
fQ3_1=figure();
plot(t_Fcalc,10^(-6)*F_drag,t_Fcalc,10^(-6)*F_inert)
xlabel('time (s)','Fontsize',14)
ylabel('Inline force [MN]','Fontsize',14)
legend({'Drag force','Inertial force'},'Fontsize',12)
title('Drag and inertial forces','Fontsize',14)
enhance_plot('TIMES',16,1.5)
saveas(fQ3_1,'Q3_1.png');

% Total inline force
fQ3_2=figure();
plot(t_Fcalc,10^(-6)*F_tot,'g',t_Fcalc,10^(-6)*F_tot_w,'k--')
xlabel('time (s)','Fontsize',14)
ylabel('Inline force [MN]','Fontsize',14)
legend({'Without wheeler stretching','With wheeler stretching'},'Fontsize',12)
title('Total inline force','Fontsize',14)
enhance_plot('TIMES',16,1.5)
saveas(fQ3_2,'Q3_2.png');

% Drag force with and without wheeler stretching
fQ3_3=figure();
plot(t_Fcalc,10^(-3)*F_drag,'g',t_Fcalc,10^(-3)*F_drag_w,'k--')
xlabel('time (s)','Fontsize',14)
ylabel('Force [kN]','Fontsize',14)
legend({'Without wheeler stretching','With wheeler stretching'},'Fontsize',12)
title('Drag force','Fontsize',14)
enhance_plot('TIMES',16,1.5)
saveas(fQ3_3,'Q3_3.png');

% Inertial force with and without wheeler stretching

fQ3_4=figure();
plot(t_Fcalc,10^(-6)*F_inert,'g',t_Fcalc,10^(-6)*F_inert_w,'k--')
xlabel('time (s)','Fontsize',14)
ylabel('Force [MN]','Fontsize',14)
legend({'Without wheeler stretching','With wheeler stretching'},'Fontsize',12)
title('Inertial force','Fontsize',14)
enhance_plot('TIMES',16,1.5)
saveas(fQ3_4,'Q3_4.png');

%% Question 4.
%according to tempel et machin with D/L=0.04 CM ~2 

%% Question 5 

% Calculate eta with the same time vector than the forces 
eta=a*cos(omega.*t_Fcalc);

% Plot the time series with the overturning moment

fQ5_1=figure();
plot(t_Fcalc,eta,t_Fcalc,10^(-6)*F_tot,t_Fcalc,10^(-6)*M/2)
xlabel('time (s)','Fontsize',14)
ylabel('','Fontsize',14)
legend({'Free surface elevation[m]','Inline force[MN]','Overturning moment [MN.m]'},'Fontsize',12)
title('Time series','Fontsize',14)
enhance_plot('TIMES',16,1.5)
saveas(fQ5_1,'Q5_1.png');

NUM = [10 15 20 25 30 35 40 50 70 90]; % trying different number of points in the vertical direction

% We will compare the amplitude of the inline force and overturning moment
% to see the effect of the number of points

F_amplitudes = zeros(1,length(NUM));
M_amplitudes = zeros(1,length(NUM));

for i =1:(length(NUM))
    [~,~,F5,M5] = FCalc(f,h,g,D,a,false,NUM(i));
    F_amplitudes(i)=max(F5);
    M_amplitudes(i)=max(M5);
end

% Plot

fQ5_2=figure();

subplot(1,2,1)
hold on
plot(NUM, 10^(-6)*F_amplitudes)
xlabel('Points along the vertical')
ylabel('F [MN]')
title('Inline force')
enhance_plot('TIMES',16,1.5)
plot(NUM,10^(-6)*F_amplitudes(end)*ones(1,length(NUM)),'k')
hold off

subplot(1,2,2)
plot(NUM, 10^(-6)*M_amplitudes)
xlabel('Points along the vertical')
ylabel('M [MN.m]')
title('Overturning moment')
enhance_plot('TIMES',16,1.5)
hold on
plot(NUM,10^(-6)*M_amplitudes(end)*ones(1,length(NUM)),'k')

saveas(fQ5_2,'Q5_2.png');


