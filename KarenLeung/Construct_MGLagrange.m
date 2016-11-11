function xcombdot =  Construct_MGLagrange(t,xcomb,Wd,WdO,rholistC, rholistO,xstar, ustar,B,C)

n = 3;                                          % Degree
xcurr = xcomb(1:n);                             % Current states
xhat = xcomb(n+1:end);                          % Estiamted states
xdiff = xstar-xhat;                             % Difference in position
Dx = (B'/Wd)*xdiff;                             % Differential
rholist = rholistC;                             % Controller rho coefficients

if (length(rholist)>1)
    const = rholist(1) + rholist(2)*xcurr(1) + rholist(3)*xcurr(2) ...
            + rholist(5)*xcurr(1)^2 + rholist(6)*xcurr(1)*xcurr(2) ...
            + rholist(7)*xcurr(2)^2 + rholist(10)*xcurr(3)^2;
    lin = rholist(2)*xdiff(1) + rholist(3)*xdiff(2) ...
          + rholist(5)*2*xcurr(1)*xdiff(1) + rholist(6)*(xcurr(2)*xdiff(1) ...
          + xcurr(1)*xdiff(2)) + rholist(7)*2*xcurr(2)*xdiff(2) ...
          + rholist(10)*2*xcurr(3)*xdiff(3);
    quad = rholist(5)*xdiff(1)^2 + rholist(6)*xdiff(1)*xdiff(2) ...
           + rholist(7)*xdiff(2)^2 + rholist(10)*xdiff(3)^2;

    rhointp = polyint([quad lin const]);
    rhoint = polyval(rhointp,1)-polyval(rhointp,0);
else
    rhoint= rholist;
end

u = ustar + 0.5*rhoint*Dx;                      % Controller

y = C*xcurr;
[k,~] = size(C);
H = [eye(n) C'
    C zeros(k)];
b = [xhat;y];

sol = H\b;
xbar = sol(1:n);
xdiff = xhat-xbar;

Dx = (WdO\C')*C*xdiff;                          % Differential
rholist = rholistO;                             % Observer rho coefficients

if (length(rholistO)>1)
    const = rholist(1) + rholist(2)*xbar(1) + rholist(3)*xbar(2) ...
            + rholist(5)*xbar(1)^2 + rholist(6)*xbar(1)*xbar(2) ...
            + rholist(7)*xbar(2)^2 + rholist(10)*xbar(3)^2;
    lin = rholist(2)*xdiff(1) + rholist(3)*xdiff(2) + rholist(5)*2*xbar(1)*xdiff(2)...
            + rholist(6)*(xbar(2)*xdiff(1) + xbar(1)*xdiff(2)) ...
            + rholist(7)*2*xbar(2)*xdiff(2) + rholist(10)*2*xbar(3)*xdiff(3);
    quad = rholist(5)*xdiff(1)^2 + rholist(6)*xdiff(1)*xdiff(2) ...
            + rholist(7)*xdiff(2)^2 + rholist(10)*xdiff(3)^2;

    rhointp = polyint([quad lin const]);
    rhoint = polyval(rhointp,1)-polyval(rhointp,0);
else
    rhoint= rholistO;
end

psiCO = 0.6;
beta = 0.5;
sigma = 7;
gamma = 0.3;

xdot = [-sigma*xcurr(1)^2 - sigma*xcurr(1)*(2*xcurr(2) + xcurr(2)^2) ; 
        -xcurr(3) - 1.5*xcurr(2)^2 - 0.5*xcurr(2)^3 - 3*xcurr(1)*xcurr(2) - 3*xcurr(1) ;
        xcurr(2) ] + B*u;



xhatdot = [-sigma*xhat(1)^2 - sigma*xhat(1)*(2*xhat(2) + xhat(2)^2) ;
           -xhat(3) - 1.5*xhat(2)^2 - 0.5*xhat(2)^3 - 3*xhat(1)*xhat(2) - 3*xhat(1) ;
           xhat(2)] + B*u- 0.5*rhoint*Dx;

xcombdot = [xdot;xhatdot];