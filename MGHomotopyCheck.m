function [Knew, kappa] = MGHomotopyCheck
syms x1s x2s lambda; 

kappa = int(integrand(lambda*x1s, lambda*x2s)*[x1s; x2s], lambda, 0, 1); 
Knew = jacobian(kappa, [x1s, x2s]); 
Knew = matlabFunction(vpa(Knew)); 
kappa = matlabFunction(vpa(kappa)); 

end
function val = integrand(x1, x2)

val = (1+x1^2 + x2^2)* [224.0878976-168.5534*x1+229.4303*x1^2+111.0209*x2^2, -64.6869246+48.6559*x1-66.2291*x1^2-32.0481*x2^2]; 

end