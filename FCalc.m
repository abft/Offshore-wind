function [F_drag,F_inert,F_tot,M]=FCalc(f,h,g,D,a,iWheeler,num)

%iWheeler being true if wheeler stretching considered, false otherwise

% Initialisation 
C_D = 1; % drag coefficient used for the monopile
C_M = 2; % inertia coefficient
rho = 1024; % sea water density
omega = 2*pi*f; 
k =kSolve(f,g,h);
A = pi*D^2/4; % cross area of the monopile

F_drag=zeros(1,200);
F_inert=zeros(1,200);
F_tot=zeros(1,200); % Inline force
M=zeros(1,200); % overturning moment

for t=1:200
    eta = a*cos(omega*t); % free surface elevation at x=0 
    z_phys = linspace(-h,eta,num);
    if iWheeler ==true
        z_calc = (z_phys - eta)/(1+eta/h);
    else
        z_calc = z_phys; % no stretching
    end
    for i=1:num
        % the horizontal velocity and acceleration are calculated with
        % z_calc
        u(i) = omega*a*cosh(k*(z_calc(i)+h))/sinh(k*h)*cos(omega*t); %airy
        dudt(i) = -omega^2*a*cosh(k*(z_calc(i)+h))/sinh(k*h)*sin(omega*t);
        dF_drag(i) = 1/2*rho*C_D*abs(u(i))*u(i);
        dF_inert(i) = rho*C_M*A*dudt(i);
        dM(i) = (dF_drag(i) + dF_inert(i))*(h + z_phys(i));
    end
    F_drag(t) = trapz(z_phys,dF_drag);
    F_inert(t) = trapz(z_phys,dF_inert); 
    M(t) = trapz(z_phys,dM);
end
F_tot = F_drag + F_inert;

end