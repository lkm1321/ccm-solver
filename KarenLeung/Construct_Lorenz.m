function xcombdot =  Construct_Lorenz(t,xcomb,Wd,WdO, rholistC, rholistO,xstar, ustar,B)
%-------------------------------------------------------------------------%
% This function finds the controller and observer design specific for the
% Coupled Lorenz Strange Attractor.                  
% Inputs are contraction metrics found, the first state estimation and the 
% coefficients of the polynomial rho found  previously. This specifically 
% works with the 3 state system because the this function is tailored with
% the constant W's and structure of rho found previously. This function 
% outputs the controlled xdots and the observed xdots.
% Author: Karen Leung
% -- Last updated 9/25/2014 -- %
%-------------------------------------------------------------------------%

n = 3;                                              % Number of state
xcurr = xcomb(1:n);                                 % Current state
xhat = xcomb(n+1:end);                              % Estimated state    
sigma = 10;
bb = 8/3;                                           % Chaotic parameters
r = 28;
C = [1 0 0;  0 0 1];                                % Observability matrix    
xdiff = xstar-xhat;                                 % Difference in position
Dx = (B'/Wd)*xdiff;                                 % Differential
rholist = rholistC;                                 % Coefficients of rho

if (length(rholist)>1)
    const = rholist(1) + rholist(5)*xcurr(1)^2 + rholist(7)*xcurr(2)^2 + rholist(9)*xcurr(2)*xcurr(3) + rholist(10)*xcurr(3)^2;
    lin = 2*(rholist(5)*xcurr(1)*xdiff(1) + rholist(7)*xcurr(2)*xdiff(2) + rholist(10)*xcurr(3)*xdiff(3)) + rholist(9)*(xcurr(2)*xdiff(3) + xcurr(3)*xdiff(2)) ;
    quad = rholist(5)*xdiff(1)^2 + rholist(7)*xdiff(2)^2 + rholist(9)*xdiff(2)*xdiff(3) + rholist(10)*xdiff(3)^2;
        
    rhointp = polyint([quad lin const]);
    rhoint = polyval(rhointp,1)-polyval(rhointp,0);
else
    rhoint= rholist;
end

u = ustar + 0.5*rhoint*Dx;                          % Controller constructed

%% Observer
y = C*xcurr;                                        % Observable states
H = [eye(n) C'
    C zeros(2)];
b = [xhat;y];
sol = H\b;
xbar = sol(1:n);
xdiff = xhat-xbar;                                  % Different in estimate    
Dx = (WdO\C')*C*xdiff;                              % Differential            
rholist = rholistO; 
if (length(rholistO)>1)

    const = rholist(1) + rholist(5)*xbar(1)^2 + rholist(7)*xbar(2)^2 + rholist(9)*xbar(2)*xbar(3) + rholist(10)*xbar(3)^2;
    lin = 2*(rholist(5)*xbar(1)*xdiff(1) + rholist(7)*xbar(2)*xdiff(2) + rholist(10)*xbar(3)*xdiff(3)) + rholist(9)*(xbar(2)*xdiff(3) + xbar(3)*xdiff(2)) ;
    quad = rholist(5)*xdiff(1)^2 + rholist(7)*xdiff(2)^2 + rholist(9)*xdiff(2)*xdiff(3) + rholist(10)*xdiff(3)^2;    

    rhointp = polyint([quad lin const]);
    rhoint = polyval(rhointp,1)-polyval(rhointp,0);
else
    rhoint= rholistO;
end

uc = 0.5*rhoint*Dx;                                  % Observer design   

% State equations
xdot = [sigma*(xcurr(2) - xcurr(1)) ; r*xcurr(1) - xcurr(3) ;  - bb* xcurr(3)] + B*u;
xhatdot = [sigma*(xhat(2) - xhat(1)) ; r*xcurr(1) - xhat(3) ;  - bb* xhat(3)] + B*u - uc;

xcombdot = [xdot;xhatdot];



