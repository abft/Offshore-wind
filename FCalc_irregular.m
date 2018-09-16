function [eta_tot,F,M]=FCalc_irregular(f,h,g,D,a,phi)

iWheeler = true;

rho = 1024;

C_D = 1;
C_M = 2;
A = pi*D^2/4;
num =40;

eta = zeros(length(f),3600);
F = zeros(1,3600);
M = zeros(1,3600);
F_drag = eta;
F_inert = eta;
T = [1:1:3600];

for j=1:length(f)
    omega = 2*pi*f(j);
    k =kSolve(f(j),g,h);
    
    for l=1:length(T)
        t=T(l);
        eta(j,l) = a(j)*cos(omega*t+phi(j));
        z_phys = linspace(-h,eta(j,l),num);
        dz = (z_phys(2)-z_phys(1))*ones(1,40); %constant
        
        if iWheeler ==true
            z_calc = (z_phys - eta(j,l))/(1+eta(j,l)/h);
        else
            z_calc = z_phys;
        end
        u = zeros(1,num);
        dudt = zeros(1,num);
        dF_drag = zeros(1,num);
        dF_inert = zeros(1,num);
        for i=1:num
            u(i) = omega*a(j)*cosh(k*(z_calc(i)+h))/sinh(k*h)*cos(omega*t+phi(j)); %airy
            dudt(i) = -omega^2*a(j)*cosh(k*(z_calc(i)+h))/sinh(k*h)*sin(omega*t+phi(j));
            dF_drag(i) = 1/2*rho*C_D*abs(u(i))*u(i);
            dF_inert(i) = rho*C_M*A*dudt(i);
            dM(i) = (dF_drag(i) + dF_inert(i))*(h + z_phys(i)); %zphys/calc
        end
        M(l) = M(l) + trapz(z_phys,dM);
        F(l) = F(l) + trapz(z_phys,dF_drag + dF_inert);
        
    end
    eta_tot = sum(eta,1);
end


end