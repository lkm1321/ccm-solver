%% Moore Greitzer 3-state system using S-Procedure, Krstic et al.
% ----------------------------------------------------------------------- %
% This function designs a outputfeedback controller for the 3-state
% Moore-Greitzer jet engine given by Krstic et al. To account for the
% restriction on the R variable, Lagrange Multipliers are used.
% Author: Karen Leung 310241847
% -- Last updated 10/6/14 -- %
% ----------------------------------------------------------------------- %
function [W,c1,WO,c2,tout,xout] = Feedback_MGLagrange()
addpath('~/mosek/7/toolbox/r2013a');

R = sdpvar(1);                      % State variable
phi = sdpvar(1);                    % State variable            
psi = sdpvar(1);                    % State variable

ep = 1e-6;
delta = sdpvar(3,1);                % Dummy variables
xi = sdpvar(3,1);                   % Dummy variables
sigma = 7;                          % Parameter

% state equation
f = [-sigma*R^2-sigma*R*(2*phi+phi^2) ;
     -psi-1.5*phi^2-0.5*phi^3-3*R*phi-3*R ;
     phi];

B = [1 1 ; 1 0 ; 0 1];              % Control matrix
C = [0 1 1 ;1 0 1];                 % Observer matrix

% Jacobian Matrix
F = jacobian(f,[R , phi , psi]);    % Jacobian matrix

%% Control Contraction
r_rho = 2;                          % Degree of rho
r_mult = 2;                         % Degree of Lagrange multipliers
lambda = 0.99;                      % Convergence rate
W = sdpvar(3,3);                    % W metric
[rho, c1, v1] = polynomial([R phi psi],r_rho,0);    % rho polynomial
rhoconstr = sos(rho);               % rho constraint

H = -(F*W + W*F'- rho*B*B'+ 2*lambda*W);    % LMI inequality
P = delta'*H*delta;                         % Quadratic for of H
lower = 1E-3;                               % W lower bound
Wconstr = [W>=lower*eye(3)];                
upper = 10;                                 % W upper bound
Wconstr2 = [W<=upper*eye(3)];

% Restriction on R bounded between 0 and 1
q1 = R;                                     % R constraint 1
q2 = 1-R;                                   % R constraint 2
% Lagrange multipliers
[lambda1, L1, v2] = polynomial([R phi psi delta'],r_mult,0);
[lambda2, L2, v3] = polynomial([R phi psi delta'],r_mult,0);
PP = P - lambda1*q1 - lambda2*q2;           % Lagrange equation
constr1 = sos(lambda1);                     % Multiplier constraint
constr2 = sos(lambda2);                     % Multiplier constraint
constr3 = sos(PP);                          % LMI constraint
% Solving SOS programming
options = sdpsettings('solver','mosek');
consList = [Wconstr, rhoconstr, constr1, constr2, constr3];
coefList=[W(:); c1; L1; L2];
% search for SOS certificate
[sol, q, Q, res] = solvesos(consList, [], options, coefList);
verifiedrho = clean(replace(rho, coefList, double(coefList)), 1e-9);
sdisplay(verifiedrho)
W = double(W)

%% Observer design
r_rho = 2;                                  % rho degree
lambda = 0.99;                              % Convergence rate
WO = sdpvar(3);                             % WO metric
[rho,c2, v2] = polynomial([R phi psi],r_rho,0); % rho polynomial
H = (F'*WO+WO*F-rho*C'*C+2*lambda*WO);      % LMI inequality
Wdotconstr = sos(-delta'*H*delta);          % LMI constraint
lower = 1E-4;                               % W lower bound
Wconstr = [WO>=lower*eye(3)];
rhoconstr = sos(rho);                       % rho constraint
% Solving SOS programming
options = sdpsettings('solver','mosek');
consList = [Wdotconstr, Wconstr, rhoconstr];
coefList=[WO(:); c2];
% search for SOS certificate
[sol, q, Q, res] = solvesos(consList, [], options, coefList);
verifiedrho = clean(replace(rho, coefList, double(coefList)), 1e-7);
sdisplay(verifiedrho)
WO = double(WO)
%% Constructing the Controller and Observer
x0 = [0.1;0.1;0.1;0.15;0.5;0.2];            % Initial conditions
xstar = [0;0;0];                            % Desired position
ustar = 0;  

rholistC = double(c1);
rholistO = double(c2);
tic
[tout,xout]= ode45(@(t,xcurr) Construct_MGLagrange(t,xcurr,W,WO,rholistC,rholistO,xstar, ustar,B,C), [0 10], x0);
toc

