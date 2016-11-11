function xcombdot =  Construct_VdPSingleControl(t,xcomb,Wd,WdO,rholistC, rholistO, ustar,mu,t_lead,x_lead)
%-------------------------------------------------------------------------%
% This function finds the controller and observer design specific for Van
% Der Pol oscillator where three oscillators need to be synchronised with
% each other using only one controller/observer
% Author: Karen Leung 310241847
% -- Last updated 10/6/14 -- %
%-------------------------------------------------------------------------%
xstar = interp1(t_lead,x_lead,t);
xstar = xstar';
n = 2;
xcurr = xcomb(1:n);
xhat = xcomb(n+1:end);

B = [0;1];
C= [1 0];

% Controller based on current state estimate xhat

xdiff = xstar-xhat;

Dx = (B'/Wd)*xdiff;

rholist = rholistC;

if (length(rholist)>1)
    const = rholist(1) + rholist(4)*xcurr(1)^2 + rholist(6)*xcurr(2)^2;
    lin = rholist(4)*2*xcurr(1)*xdiff(1) + rholist(6)*2*xcurr(2)*xdiff(2);
    quad = rholist(4)*xdiff(1)^2 + rholist(6)*xdiff(2)^2;

    rhointp = polyint([quad lin const]);
    rhoint = polyval(rhointp,1)-polyval(rhointp,0);
else
    rhoint= rholist;
end

u = ustar + 0.5*rhoint*Dx;

% Observer to adjust xhat based on y


y = C*xcurr;

H = [eye(n) C'
    C 0];
b = [xhat;y];

sol = H\b;
xbar = sol(1:n);

xdiff = xhat-xbar;

Dx = (WdO\C')*C*xdiff;

rholist = rholistO;

if (length(rholistO)>1)
    const = rholist(1) + rholist(4)*xbar(1)^2 +rholist(6)*xbar(2)^2 ...
            + rholist(13)*xbar(1)^2*xbar(2)^2 + rholist(11)*xbar(1)^4 ...
            + rholist(15)*xbar(2)^4;
    lin = rholist(4)*2*xbar(1)*xdiff(1) + rholist(6)*xbar(2)*xdiff(2) ...
          + rholist(13)*2*(xbar(2)*xdiff(2)*xbar(1)^2 + xbar(1)*xdiff(1)*xbar(2)^2) ...
          + rholist(11)*4*xbar(1)^2*xdiff(1) +rholist(15)*xbar(2)^3*xdiff(2);
    quad = rholist(4)*xdiff(1)^2 + rholist(6)*xdiff(2)^2 ... 
          + rholist(13)*(4*xbar(1)*xbar(2)*xdiff(1)*xdiff(2) ...
          + xbar(2)^2*xdiff(1)^2 + xbar(1)^2*xdiff(2)^2) ...
          + rholist(11)*6*xbar(1)^2*xdiff(1)^2 + rholist(15)*6*xbar(2)^2*xdiff(2)^2;
    cube = rholist(13)*2*(xbar(1)*xdiff(1)*xdiff(2)^2 + xbar(2)*xdiff(2)*xdiff(1)^2) ...
          + rholist(11)*4*xbar(1)*xdiff(1)^3 + rholist(15)*4*xbar(2)*xdiff(2)^3;
    quar = rholist(13)*xdiff(1)^2*xdiff(2)^2 + rholist(11)*xdiff(1)^4 ...
          + rholist(15)*xdiff(2)^4;

    rhointp = polyint([quar cube quad lin const]);
    rhoint = polyval(rhointp,1)-polyval(rhointp,0);
else
    rhoint= rholistO;
end

uc = 0.5*rhoint*Dx; 

    xdot =[mu*(xcurr(2) - (xcurr(1)^3/3 - xcurr(1))); -xcurr(1)/mu + u] ; % Where you actually are
    xhatdot = [mu*(xhat(2) - (xhat(1)^3/3 - xhat(1))) ; -xhat(1)/mu + u] - uc; % Where you observe to be


xcombdot = [xdot;xhatdot];