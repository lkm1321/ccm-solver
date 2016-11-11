function MGIntegrableCheck
clc;
clear; 
close all; 

[~, kappa] = MGHomotopyCheck; 

MGStabilize = @(t, x)(MGmodel(t, x, 0, kappa)); 

[t, x] = ode45(MGStabilize, [0 10], [20 20]); 
figure, plot(t, x(:,1));
figure, plot(t, x(:,2)); 

x1gv = linspace(-20, 20, 100); 
x2gv = linspace(-20, 20, 100); 
[x1g, x2g] = meshgrid(x1gv, x2gv); 
z = NaN(100, 100); 

for i = 1:100
    for j = 1:100
       z(i, j) = kappa(x1gv(i), x2gv(j));  
    end
end
figure, surf(x1g, x2g, z); xlabel('x1'), ylabel('x2'), zlabel('u'); 
end

function xdot = MGmodel(~, xcurr, xstar1, kappa)

    x1 = xcurr(1);
    x2 = xcurr(2); 

    % For every xstar1, there is a stable xstar2. 
    xstar2 =  -0.5*xstar1^3 - 1.5*xstar1^2;
    
    u = kappa(x1, x2); 
       
    ulower = kappa(xstar1, xstar2);  
        
    ustar = -xstar1; 
         


    xdot =[-0.5*x1^3-1.5*x1^2-x2;... 
                x1 + u + ulower + ustar]; 

end 