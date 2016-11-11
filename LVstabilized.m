%% This controller stabilizes to manifold 

function out = LVstabilized(t, x)

x1 = x(1);
x2 = x(2); 

x1star = x(3); 
x2star = x(4); 

%e1 = x1 - x1star;
%e2 = x2 - x2star; 
% u = - (5.2425*x1-0.0065*x2-1.0185*x1*x2+0.7484*x1^2+0.1678*x2^2-0.2100*x1^2*x2+0.7093*x1*x2^2+2.3617*x1^3-2.4886e-04*x2^3); 
% u = - (44.9856*x1+6.0912*x2-0.8021*x1*x2-0.4687*x1^2); 
u = - (44.9856*x1+6.0912*x2-0.8021*x1*x2-0.4687*x1^2); 
ulower =  + (44.9856*x1star+6.0912*x2star-0.8021*x1star*x2star-0.4687*x1star^2); 
%  u = - (44.9856*e1+6.0912*e2-0.8021*e1*e2-0.4687*e1^2); 
% u = - (5.2425*x1-0.0065*x2-1.0185*x1*x2+0.7484*x1^2+0.1678*x2^2-0.2100*x1^2*x2+0.7093*x1*x2^2+2.3617*x1^3-2.4886e-04*x2^3); 
% u = 0; 

xdot = [x1-x1*x2 + u + ulower;
            -x2+x1*x2]; 
xstardot = [x1star-x1star*x2star;
                   -x2star+x1star*x2star];
out = [xdot(:); xstardot(:)]; 
               
end 