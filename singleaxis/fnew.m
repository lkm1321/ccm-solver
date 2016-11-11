function xdot = fnew(x,mu)

global k;

if (isempty(k)) 
   k = 500*10^-6/0.005;  
end

theta = x(1);
omega = x(2); 

%mu = controller(theta,omega,k); 

% disp(theta);

xdot = zeros(2,1);
xdot(1) = omega;
xdot(2) = mu*sin(theta)*k; 


end
