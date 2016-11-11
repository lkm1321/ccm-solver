function xcombdot =  Construct_MG(t,xcomb,Wd,WdO,rholistC, rholistO,xstar, ustar)
%-------------------------------------------------------------------------%
% This function finds the controller and observer design specific for the
% 2-state Moore-Greitzer model. Inputs are contraction metrics found, the
% first state estimation and the coefficients of the polynomial rho found
% previously. This specifically works with the 2 state system because the
% this function is tailored with the constant W's and structure of rho
% found previously. This function outputs the controlled xdots and the
% observed xdots.
% Author: Ian Manchester
% -- Last updated 1/5/2014 -- %
%-------------------------------------------------------------------------%

n = 2;                                      % Number of states
xcurr = xcomb(1:n);                         % Current state
xhat = xcomb(n+1:end);                      % Estimated state    
B = [0;1];                                  % Controlability matrix
C = [0 1];                                  % Observability matrix

xdiff = xstar-xhat;                         % Difference in position
Dx = (B'/Wd)*xdiff;                         % Differential component
rholist = rholistC;                         % Coefficient of rho controller

if (length(rholist)>1)
    const = rholist(1)+rholist(2)*xcurr(1)+rholist(4)*xcurr(1)^2+rholist(6)*xcurr(2)^2;
    lin = rholist(2)*xdiff(1)+2*rholist(4)*xdiff(1)*xcurr(1)+2*rholist(6)*xdiff(2)*xcurr(2);
    quad = rholist(4)*xdiff(1)^2+rholist(6)*xdiff(2)^2;

    rhointp = polyint([quad lin const]);
    rhoint = polyval(rhointp,1)-polyval(rhointp,0);
else
    rhoint= rholist;
end

u = ustar + 0.5*rhoint*Dx;                  % State feedback controller

y = C*xcurr;                                % Obsevability equation
H = [eye(n) C'
    C 0];
b = [xhat;y];
sol = H\b;
xbar = sol(1:n);
xdiff = xhat-xbar;                          % Difference in position
Dx = (WdO\C')*C*xdiff;                      % Differential
rholist = rholistO;

if (length(rholistO)>1)
    const = rholist(1)+rholist(2)*xbar(1)+rholist(4)*xbar(1)^2+rholist(6)*xbar(2)^2;
    lin = rholist(2)*xdiff(1)+2*rholist(4)*xdiff(1)*xbar(1)+2*rholist(6)*xdiff(2)*xbar(2);
    quad = rholist(4)*xdiff(1)^2+rholist(6)*xdiff(2)^2;

    rhointp = polyint([quad lin const]);
    rhoint = polyval(rhointp,1)-polyval(rhointp,0);
else
    rhoint= rholistO;
end

uc = 0.5*rhoint*Dx;                         % Observer

% State equation of current state
xdot =[-0.5*xcurr(1)^3-1.5*xcurr(1)^2-xcurr(2);xcurr(1)+u ]; 
% State equation of estimated state
xhatdot = [-0.5*xhat(1)^3-1.5*xhat(1)^2-xhat(2);xhat(1)+u ]-uc;
xcombdot = [xdot;xhatdot];