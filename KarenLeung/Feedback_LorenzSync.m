%% Contraction Analysis via SOS programming Partial Sync of Lorenz System
%-------------------------------------------------------------------------%
% This code finds a controller for the Lorenz Equations. Its parameters are
% the critical case where the system is chaotic. It uses Partial
% Contraction as the xstar input is the trajectory of another chaotic
% system that starts in a slightly different position.
% Author: Karen Leung 310241847
% -- Last updated 10/6/2014 -- %
%-------------------------------------------------------------------------%
function [W,c1,WO,c2,tout,xout,tlead,xlead,tback,xback] = Feedback_LorenzSync()
addpath 'c:\Program Files\mosek\7\toolbox\r2013a'

%% Setting up variables and parameters
sdpvar x y z                                    % State variables
sigma = 10;                                     % Chaotic parameter
b = 8/3;                                        % Chaotic parameter
r = 28;                                         % Chaotic parameter
delta = sdpvar(3,1);                            % Dummy variable
% State equation
xdot = sigma*(y - x);
ydot =  r*x - y - x*z;
zdot = x*y - b*z;

B = [1 0 ;0 1 ;0 0];                            % Control matrix
C = [1 0 0;  0 0 1];                            % Observer matrix    
f = [xdot ; ydot ; zdot];                       % Setting up the EOM
F = jacobian(f , [x , y , z]);                  % Calculating Jacobian

%% Contraction Metric set up
W = sdpvar(3,3);                                % W metric
r_rho = 2;                                      % rho degree
lambda = 2.66;                                  % Convergence rate
[rho,c1, v1] = polynomial([x y z],r_rho,0);     % rho polynomial
rhoconstr = sos(rho);                           % rho constraint
H = (F*W + W*F' - rho*B*B' + 2*lambda*W);       % LMI inequality             
Wdotconstr = sos(-delta'*H*delta);              % LMI constraint

% Bounds
lower = 0.1;
Wconstr = [W>=lower*eye(3)];
upper = 5;
Wconstr2 = [W<=upper*eye(3)];

% Solving
options = sdpsettings('solver','mosek');
consList = [Wdotconstr, Wconstr, Wconstr2,rhoconstr];
coefList=[c1;W(:)];
% search for SOS certificate
[sol, q, Q, res] = solvesos(consList, [], options, coefList);
verifiedrho = clean(replace(rho, coefList, double(coefList)), 1e-7);
sdisplay(verifiedrho)
W = double(W)

%% Observer Contraction Metric
r_rho = 2;                                      % rho degree
lambda = 47.6;                                  % Convergence rate
WO = sdpvar(3,3);                               % WO metric
[rho,c2, v2] = polynomial([x y z],r_rho,0);     % rho polynomial
H = (F'*WO + WO*F - rho*C'*C + 2*lambda*WO);    % LMI inequality
Wdotconstr = sos(-delta'*H*delta);              % LMI constraint
lower = 0.01;
Wconstr = [WO>=lower*eye(3)];                   % W lower bound
rhoconstr = sos(rho);                           % rho constraint
% SOS programming solving
options = sdpsettings('solver','mosek');
consList = [Wdotconstr, Wconstr, rhoconstr];
coefList=[WO(:); c2];
% search for SOS certificate
[sol, q, Q, res] = solvesos(consList, [], options, coefList);
verifiedrho = clean(replace(rho, coefList, double(coefList)), 1e-7);
sdisplay(verifiedrho)
WO = double(WO)


%% Constructing Controller

x0 = [5 ; 2 ; 3 ; 5.4 ; 2.5 ; 2.5];

% xstar = [0;0;0];          % where you want to be  
ustar = 0;               

rholistC = double(c1);
rholistO = double(c2);
t_max = 5;
dpar = 1.5;
lead = @(t,x) [sigma*dpar*(x(2) - x(1)); r*dpar*x(1) - x(2) - x(1)*x(3);x(1)*x(2) - b*dpar*x(3)];
[tlead,xlead] = ode45(lead,[0,t_max],x0(1:3)*1.7);
x_noise = .05.*randn(length(xlead),3).*xlead + xlead ;
t_on = 2;
back = @(t,x) [sigma*(x(2) - x(1)); r*x(1) - x(2) - x(1)*x(3);x(1)*x(2) - b*x(3)];
[tback,xback] = ode45(back,[0,t_on],x0(1:3));
x_on = xback(end,:);

tic
[tout,xout]= ode45(@(t,xcurr) Construct_LorenzSync(t,xcurr,W,WO,rholistC,rholistO,ustar,B,tlead,xlead), [t_on t_max], [x_on';x_on'*1.5]);
toc







