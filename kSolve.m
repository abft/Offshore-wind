function k =kSolve(f,g,h)

omega = 2*pi.*f;
k0 = sqrt(omega.^2/g);
k = fzero( @(k) omega.^2 -g*k*tanh(k*h),k0);

end