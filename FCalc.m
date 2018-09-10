function [F_drag,F_inert,F_tot,M]=FCalc(f,h,g,D,a,iWheeler,num)
% f = 0.0833;
% h =30;
% g = 9.81;
% D = 6;
% a = 3;
% iWheeler = true;

%iWheeler being true if wheeler stretching considered, false otherwise
rho = 1024;
omega = 2*pi*f;
k =kSolve(f,g,h);
C_D = 1;
C_M = 2;
A = pi*D^2/4;

for t=1:200
    eta = a*cos(omega*t);
    z_phys = linspace(-h,eta,num);
    dz = (z_phys(2)-z_phys(1))*ones(1,40); %constant
    if iWheeler ==true
        z_calc = (z_phys - eta)/(1+eta/h);
    else
        z_calc = z_phys;
    end
    for i=1:num
        u(i) = omega*a*cosh(k*(z_calc(i)+h))/sinh(k*h)*cos(omega*t); %airy
        dudt(i) = -omega^2*a*cosh(k*(z_calc(i)+h))/sinh(k*h)*sin(omega*t);
        dF_drag(i) = 1/2*rho*C_D*abs(u(i))*u(i);
        dF_inert(i) = rho*C_M*A*dudt(i);
        dM(i) = (dF_drag(i) + dF_inert(i))*(h + z_calc(i));
    end
    F_drag(t) = trapz(z_phys,dF_drag);
    F_inert(t) = trapz(z_phys,dF_inert); 
    M(t) = trapz(z_phys,dM);
end
F_tot = F_drag + F_inert;


end