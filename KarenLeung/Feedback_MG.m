%% Moore Greitzer 2-state system, Krstic et al.
%-------------------------------------------------------------------------%
% This code uses the model of the Moore-Greitzer engine provided by Krstic
% et al. Using the theory described by Manchester and Slotine in 'Control
% Contraction Metrics and Universal Stabilizability', it can be proven that
% this system is universally stabilisable via constraction theory and a 
% state feedback controller is found to stabilise the system.
% Author: Ian Manchester
% -- Last updated 1/5/2014 -- %
%-------------------------------------------------------------------------%
function [W,c1,WO,c2,tout,xout] = Feedback_MG()

addpath '/home/lkm1321/mosek/7/toolbox/r2013a'; 
addpath(genpath('/home/lkm1321/yalmip')); 

x1 = sdpvar(1);                             % Defining variables
x2 = sdpvar(1);                             % Defining variables
ep = 1e-6;                                  % Defining tolerences
delta = sdpvar(2,1);                        % Dummy variable

x = [x1 x2]';                               % State space

f = [-0.5*x1^3-1.5*x1^2-x2;x1];             % Dynamical system
B = [0;1];                                  % Controlablity matrix
C = [0 1];                                  % Observibility matrix

F = jacobian(f,[x1 x2]);                    % Jacobian matrix


%% Contraction Metric set up
% degree of polynomial in W
r_rho = 2;                                  % Degree for rho
lambda = 1.961;                             % Convergence rate
W = sdpvar(2,2);                            % Metric
[rho,c1, v1]=polynomial([x1 x2],r_rho,0);   % rho polynomial
rhoconstr = sos(rho);                       % rho constraint
H=(F*W+W*F'-rho*B*B'+2*lambda*W);           % LMI inequality
Wdotconstr = sos(-delta'*H*delta);          % LMI constraint    
lower = 0.1;                                % Lower bound on W
Wconstr = [W>=lower*eye(2)];                
upper = 5;                                  % Upper bound on W
Wconstr2 = [W<=upper*eye(2)];

options = sdpsettings('solver','mosek');
consList = [Wdotconstr, Wconstr, Wconstr2,rhoconstr];
coefList=[c1 ; W(:)];
% search for SOS certificate
[sol, q, Q, res] = solvesos(consList, [], options, coefList);
verifiedrho = clean(replace(rho, coefList, double(coefList)), 1e-7);
sdisplay(verifiedrho); 
Wd = double(W);

%% Observer design
r_rho = 2;                                  % Degree for rho
lambda = 9.624;                             % Convergence rate
WO = sdpvar(2,2);                           % Metric
[rho,c2, v2] = polynomial([x1 x2],r_rho,0); % rho polynomial
rhoconstr = sos(rho);                       % rho constraint
H=(F'*WO+WO*F-rho*C'*C+2*lambda*WO);        % LMI inequality
Wdotconstr = sos(-delta'*H*delta);          % LMI constraint
lower = 0.01;                               % Lower bound on W
Wconstr = [WO>=lower*eye(2)];

options = sdpsettings('solver','mosek');
consList = [Wdotconstr, Wconstr, rhoconstr];
coefList=[c2 ; WO(:)];
% search for SOS certificate
[sol, q, Q, res] = solvesos(consList, [], options, coefList);
verifiedrho = clean(replace(rho, coefList, double(coefList)), 1e-7);
sdisplay(verifiedrho)
WdO = double(WO)

%% Controller and Observer construction
x0 = [-2.5;-0.1;-1.5;-0.5];                 % Current and estimated states
xstar = [0;0];                              % Desired state   
ustar = 0;                                  % Nominal control input

rholistC = double(c1);
rholistO = double(c2);
tic
[tout,xout]= ode45(@(t,xcurr) Construct_MG(t,xcurr,Wd,WdO,rholistC,rholistO,xstar, ustar), [0 10], x0);
toc
