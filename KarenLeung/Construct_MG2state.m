function xcombdot =  Construct_MG2state(t,xcomb,Wd,WdO,rholistC, rholistO,xstar, ustar)
%-------------------------------------------------------------------------%
% This function finds the controller and observer design specific for the
% 3-state Moore-Greitzer model using the 2-state approximation.
% Inputs are contraction metrics found, the first state estimation and the
% coefficients of the polynomial rho found previously. 
% This specifically works with the 2 state system because the
% this function is tailored with the constant W's and structure of rho
% found previously. This function outputs the controlled xdots and the
% observed xdots.
% Author: Karen Leung
% -- Last updated 9/3/2014 -- %
%-------------------------------------------------------------------------%

sigma = 7;                            % Parameter used for f(x)
xcurr = xcomb(1:3);                   % Takes in current state   
xhat = xcomb(4:end);                  % Takes in estimated state 
B = [0;1];                            % B matrix for nominal model  
C = [0 1];                            % C matrix for nominal model  
%% Controller based on current 2-state estimate - xhat
xdiff = xstar-xhat;                   % Difference from where you are and
                                      % where you think you are                               
Dx = (B'/Wd)*xdiff;                   % Differential - part of formula
rholist = rholistC;                   % Uses coefficients of rho 
% This separates into constant, linear and quadrative coefficients in s
if (length(rholist) > 1)
    const = rholist(1) + rholist(2)*xcurr(1) + rholist(4)*xcurr(1)^2 ...
            + rholist(6)*xcurr(2)^2;
    lin = rholist(2)*xdiff(1) + 2*rholist(4)*xdiff(1)*xcurr(1) ...
            + 2*rholist(6)*xdiff(2)*xcurr(2);
    quad = rholist(4)*xdiff(1)^2 + rholist(6)*xdiff(2)^2;

    rhointp = polyint([quad lin const]);             % construct polynomial
    rhoint = polyval(rhointp,1)-polyval(rhointp,0);  % Integrates
else
    rhoint = rholist;                                % Else constant rho
end
u = ustar + 0.5*rhoint*Dx;                     % Constructs the controller 

%% Observer to adjust 2-state xhat based on y
y = C*xcurr(1:2);                   % Takes in the first two variables only

H = [eye(2) C';
     C 0];
b = [xhat;y];
sol = H\b;
xbar = sol(1:2);
xdiff = xhat-xbar;                  % Using formula from paper
Dx = (WdO\C')*C*xdiff;              % Construct the differential
rholist = rholistO;                 % Takes in coefficients of rhoO
% This separates into constant, linear and quadrative coefficients in s
if (length(rholistO) > 1)
    const = rholist(1) + rholist(2)*xbar(1) + rholist(4)*xbar(1)^2 ...
             + rholist(6)*xbar(2)^2;
    lin = rholist(2)*xdiff(1) + 2*rholist(4)*xdiff(1)*xbar(1) ...
            + 2*rholist(6)*xdiff(2)*xbar(2);
    quad = rholist(4)*xdiff(1)^2 + rholist(6)*xdiff(2)^2;

    rhointp = polyint([quad lin const]);
    rhoint = polyval(rhointp,1)-polyval(rhointp,0);
else
    rhoint= rholistO;               % Else rhoO is constant
end
uc = 0.5*rhoint*Dx;                 % Construct the observer
% Using the nominal controller and observer, plug it back into the 3-state 
% equation and see how all 3 states behave.
xdot = [ -0.5*xcurr(1)^3 - 1.5*xcurr(1)^2 - xcurr(2) ...
            - 3*xcurr(3)*xcurr(1) - 3*xcurr(3);
        xcurr(1)+u;
        -sigma*xcurr(3)^2-sigma*xcurr(3)*(2*xcurr(1) + xcurr(1)^2)];     
xhatdot = [-0.5*xhat(1)^3-1.5*xhat(1)^2-xhat(2);xhat(1)+u ]-uc; 

xcombdot = [xdot;xhatdot];