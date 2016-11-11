clear;
clc;
close all; 

%% Add MOSEK and YALMIP Path : Please modify appropriately. 

addpath '/home/lkm1321/mosek/7/toolbox/r2013a'; 
addpath(genpath('/home/lkm1321/yalmip')); 

%% Finding CCM. Taken from Feedback_MG.m in Karen Leung's thesis

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
W = sdpvar(2,2, 'symmetric');                            % Metric. Assumed constant. 
[rho,c1, v1]=polynomial([x1 x2],r_rho,0);   % rho polynomial
rhoconstr = sos(rho);                       % rho constraint

% The below two lines are too hard for MOSEK. 
% integrability = sos(delta'*B*jacobian(rho, [x1, x2])*W*delta); 
% integrability = norm( jacobian( rho * B' / double(W), [x1, x2] ) - jacobian( rho * B' / double(W), [x1, x2])' ); 

H=(F*W+W*F'-rho*(B*B')+2*lambda*W);           % LMI inequality

% Solving for W and K at the same time is also too hard for MOSEK. 
% H = (F*W + W*F' - B * K * W - W * K' * B' + 2 * lambda * W); 

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

% This is the CCM-based differential controller. 
K = -0.5*verifiedrho * B' / double(W); 
jacK = jacobian(K', [x1, x2]); 

%% D-K type iteration for integrability certificate
% No iteration was needed in fact. Yet, to be fair, MG model was differentially flat.  
% Interesting question is, if du = Kdx ensures contraction, can we have du = (K +
% Delta) *dx a stabilizing controller for bounded Delta?  

% Fix W. 
Wd = double(W);
% The non-differential controller 
[K_integral, integralCoeffList, K_monomials] = polynomial([x1, x2], max(degree(K, [x1,x2])) + 2, 1); 
% Differential controller. Enforce integrability. 
K_differential = jacobian(K_integral, [x1,x2]);
% Maximize lambda. Resulting value is quite similar to the initial lambda
% (initial: 1.961, D-K: 1.9642)
lambda_free = sdpvar(1); 
% W_dot constraint. 
H_K = (F*Wd + Wd*F' - B*K_differential*Wd - Wd*K_differential'*B' + 2*lambda_free*Wd ); 
H_K_constr = sos(-delta'*H_K*delta); 
% Thought this might be healthy. 
lambdaconstr = [lambda_free >= (lambda-1)];  
% Constraints list. 
constList_K = [H_K_constr, lambdaconstr]; 
% Parmaeters list
coefList_free = [integralCoeffList(:); lambda_free];
% Minimize - lambda_free (minus lambda_free), subject to contraction
% condition. Allowed to vary integral coefficients, and lambda_free. 
[sol_free, q_free, Q_free, res] = solvesos(constList_K, -lambda_free, options, coefList_free); 
% Since K_integral and K_differential share coefficients list, below clean
% statements give the integrated controller. integrableK is the
% differential version, and integralK is the integrated version. 
integrableK = clean(replace(K_differential, coefList_free, double(coefList_free)), 1e-7); 
integralK = clean(replace(K_integral, coefList_free, double(coefList_free)), 1e-7); 
sdisplay(integrableK); 
sdisplay(integralK); 
% For convenience, I copy and pasted result of sdisplay(integralK) in
% MGmodel. The last argument is ustar. Because of nonlinearity, for every
% ustar there is associated equilibrium. 