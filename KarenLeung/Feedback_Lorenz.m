%% Contraction Analysis via SOS programming Moore Greitzer Model
%-------------------------------------------------------------------------%
% This code finds a controller for the Lorenz Equations. Its parameters are
% the critical case where the system is chaotic. 
% Author: Karen Leung 310241847
% -- Last updated 9/25/2014 -- %
%-------------------------------------------------------------------------%
function [W,c1,WO,c2,tout,xout] = Feedback_Lorenz(B)
addpath 'c:\Program Files\mosek\7\toolbox\r2013a'

%% Setting up variables and parameters
sdpvar x y z                                        % State variables
sigma = 10;                                         % Chaotic parameter
b = 8/3;                                            % Chaotic parameter    
r = 28;                                             % Chaotic parameter
delta = sdpvar(3,1);                                % Dummy variable

%% The system equation and its jacobian
xdot = sigma*(y - x);
ydot =  r*x - y;                                    % Equations of motion
zdot =  - b*z;
C = [1 0 0;  0 0 1];                                % Observability matrix
f = [xdot ; ydot ; zdot];                           % Setting up the EOM
F = jacobian(f , [x , y , z]);                      % Jacobian matrix
%% Controller design
W = sdpvar(3,3);                                    % Metric
r_rho = 0;                                          % Degree of rho
lambda = 2.66;                                      % Convergence rate
[rho,c1, v1] = polynomial([x y z],r_rho,0);         % rho polynomial
rhoconstr = sos(rho);                               % rho constraint
H=(F*W+W*F'-rho*B*B'+2*lambda*W);                   % LMI inequality
Wdotconstr = sos(-delta'*H*delta);                  % LMI constraint
lower = 0.1;                                        % Lower bound on W
Wconstr = [W>=lower*eye(3)];
upper = 5;                                          % Upper bound on W
Wconstr2 = [W<=upper*eye(3)];

options = sdpsettings('solver','mosek');
consList = [Wdotconstr, Wconstr, Wconstr2,rhoconstr];
coefList=[c1;W(:)];
% search for SOS certificate
[sol, q, Q, res] = solvesos(consList, [], options, coefList);
verifiedrho = clean(replace(rho, coefList, double(coefList)), 1e-7);
sdisplay(verifiedrho)
Wd = double(W)

%% Observer design
r_rho = 0;                                          % Degree of rho
lambda = 75.5;                                      % Convergence rate
WO = sdpvar(3,3);                                   % Metric
[rho,c2, v2] = polynomial([x y z],r_rho,0);         % rho polynomial
rhoconstr = sos(rho);                               % rho constraint
H=(F'*WO+WO*F-rho*C'*C+2*lambda*WO);                % LMI inequality
Wdotconstr = sos(-delta'*H*delta);                  % LMI constraint
lower = 0.01;                                       % Lower bound on WO
Wconstr = [WO>=lower*eye(3)];
                              
options = sdpsettings('solver','mosek');
consList = [Wdotconstr, Wconstr, rhoconstr];
coefList=[WO(:); c2];
% search for SOS certificate
[sol, q, Q, res] = solvesos(consList, [], options, coefList);
verifiedrho = clean(replace(rho, coefList, double(coefList)), 1e-7);
sdisplay(verifiedrho)
WdO = double(WO)

%% Constructing Controller
x0 = [5 ; 7 ; 3 ; 5.4 ; 7.5 ; 2.5];                 % Current and estimated state
xstar = [0;0;0];                                    % Desired state    
ustar = 0;                                          % Nominal control
rholistC = double(c1);                              % rho controller coefficient    
rholistO = double(c2);                              % rho observer coefficient
tic
[tout,xout]= ode45(@(t,xcurr) Construct_Lorenz(t,xcurr,Wd,WdO,rholistC,rholistO,xstar,ustar,B), [0 5], x0);
