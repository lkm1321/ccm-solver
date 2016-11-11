%% Moore Greitzer 3-state system which is turned into 2 state, Krstic et al.
%-------------------------------------------------------------------------%
% This code models the 3 state system by using the 2 state system, but
% modifying it to account for the third state. (It was infeasbile by
% directly applying the third state to the system.)  This approach is
% robust control with an uncertainty parameter involved.
% Author: Karen Leung 310241847
% -- Last updated 9/25/2014 -- %
%-------------------------------------------------------------------------%
function [W,c1,WO,c2,tout,xout] = Feedback_MG3D()
addpath 'c:\Program Files\mosek\7\toolbox\r2013a'

x1 = sdpvar(1);                                 % Defining state variable
x2 = sdpvar(1);                                 % Defining state variable
ep = 1e-6;                                      % Tolerences
delta = sdpvar(2,1);                            % Dummy variable
x = [x1 x2]';                                   % state space

f = [-0.5*x1^3-1.5*x1^2-x2;x1];                 % System of equations
L = [-3 , 0 ; 0 , 0 ];                          % Uncertainty parameter
B = [0;1];                                      % Controlability matrix
C = [0 1];                                      % Observability matrix

F1 = jacobian(f,[x1 x2]);                       % R=0 Jacobian        
F2 = jacobian(f,[x1 x2]) + L;                   % R=1 Jacobian

%% Contraction Metric set up
r_rho = 2;                                      % Degree of rho
lambda = 1.962;                                 % Convergence rate
W = sdpvar(2,2);                                % Metric    
[rho , c1 , v1] = polynomial([x1 x2],r_rho,0);  % rho polynomial
rhoconstr = sos(rho);                           % rho constraint
H1 = (F1*W+W*F1'-rho*B*B'+2*lambda*W);          % LMI inequality
H2 = (F2*W+W*F2'-rho*B*B'+2*lambda*W);          % LMI inequality
Wdotconstr1 = sos(-delta'*H1*delta);            % LMI constraint
Wdotconstr2 = sos(-delta'*H2*delta);            % LMI constraint
lower = 0.1;                                    % Lower bound on W
Wconstr = [W>=lower*eye(2)];                    
upper = 5;                                      % Upper bound on W
Wconstr2 = [W<=upper*eye(2)];

options = sdpsettings('solver','mosek');
consList = [Wdotconstr1, Wdotconstr2, Wconstr, Wconstr2,rhoconstr];
coefList=[W(:); c1];
% search for SOS certificate
[sol, q, Q, res] = solvesos(consList, [], options, coefList);
verifiedrho = clean(replace(rho, coefList, double(coefList)), 1e-7);
sdisplay(verifiedrho)
Wd = double(W)

%% Observer design
r_rho = 2;                                      % Degree of rho
lambda = 8.71;                                  % Convergence rate
WO = sdpvar(2,2);                               % Metric
[rho,c2, v2]=polynomial([x1 x2],r_rho,0);       % rho polynomial
rhoconstr = sos(rho);                           % rho constraint
H1 = (F1'*WO+WO*F1-rho*C'*C+2*lambda*WO);       % LMI inequality
H2 = (F2'*WO+WO*F2-rho*C'*C+2*lambda*WO);       % LMI inequality
Wdotconstr1 = sos(-delta'*H1*delta);            % LMI constraint
Wdotconstr2 = sos(-delta'*H2*delta);            % LMI constraint
lower = 0.01;                                   % Lower bound on WO
Wconstr = [WO>=lower*eye(2)];

options = sdpsettings('solver','mosek');
consList = [Wdotconstr1, Wdotconstr2, Wconstr, rhoconstr];
coefList=[WO(:); c2];
% search for SOS certificate
[sol, q, Q, res] = solvesos(consList, [], options, coefList);
verifiedrho = clean(replace(rho, coefList, double(coefList)), 1e-7);
sdisplay(verifiedrho)
WdO = double(WO)

%% Controller and Observer construction
x0 = [-3;-0.4;0.4;-2.5;-0.5];                   % Current and estimated state
xstar = [0;0];                                  % Desired state
ustar = 0;                                      % Nominal control
rholistC = double(c1);                          % rho control coefficients
rholistO = double(c2);                          % rho observer coefficients
tic
[tout,xout]= ode45(@(t,xcurr) Construct_MG2state(t,xcurr,Wd,WdO,rholistC,rholistO,xstar, ustar), [0 10], x0);
toc
