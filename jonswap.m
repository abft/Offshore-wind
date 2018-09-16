function [f,a,S]=jonswap(Hs, Tp, df, fHighCut, gammaJS)

    f=[df:df:fHighCut];
    fp = 1/Tp;
    r = f./fp; %ratio
    S = zeros(1,length(f));
    
    for i =1:length(f)
        if f(i) <= fp
            sigma = 0.07;
        else
            sigma = 0.09;
        end
        S(i) = 0.3125*Hs^2*Tp*r(i)^(-5)* exp(-1.25*r(i)^(-4))*(1-0.287*log(gammaJS)) ...
            *gammaJS^(exp(-0.5*((r(i)-1)/sigma)^2));
    end
    a = sqrt(2.*S.*df);
    
end
