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
lambdaO = 9.624; 
W = sdpvar(2,2, 'symmetric');                            % Metric. Assumed constant. 
M = sdpvar(2,2, 'symmetric'); 

[rho,c1, v1]=polynomial([x1 x2],r_rho,0);   % rho polynomial
rhoconstr = sos(rho);                       % rho constraint

[rhoObs, cObs, ~] = polynomial([x1, x2], r_rho, 0); 
rhoObsConstr = sos(rhoObs); 

% The below two lines are too hard for MOSEK. 
% integrability = sos(delta'*B*jacobian(rho, [x1, x2])*W*delta); 
% integrability = norm( jacobian( rho * B' / double(W), [x1, x2] ) - jacobian( rho * B' / double(W), [x1, x2])' ); 

H=(F*W+W*F'-rho*(B*B')+2*lambda*W);           % LMI inequality
Ho = (M*F+F'*M-rhoObs*(C'*C) + 2*lambdaO*M); 

% Solving for W and K at the same time is also too hard for MOSEK. 
% H = (F*W + W*F' - B * K * W - W * K' * B' + 2 * lambda * W); 

Wdotconstr = sos(-delta'*H*delta);          % LMI constraint    
lower = 0.1;                                % Lower bound on W
Wconstr = [W>=lower*eye(2)];                
upper = 5;                                  % Upper bound on W
Wconstr2 = [W<=upper*eye(2)];
options = sdpsettings('solver','mosek', 'verbose', 2);
consList = [Wdotconstr, Wconstr, rhoconstr];
coefList=[c1 ; W(:)];
% search for SOS certificate
[sol, q, Qs, res] = solvesos(consList, [], options, coefList);
verifiedrho = clean(replace(rho, coefList, double(coefList)), 1e-7);

Mdotconstr = sos(-delta'*(Ho)*delta);
Mlower = 0.01; 
Mupper = 100;
Mconstr = [ (M >= Mlower*eye(2) ), ( M <= Mupper*eye(2) )]; 

consList = [Mdotconstr, Mconstr, rhoObsConstr]; 
obsCoeffList = [cObs; M(:)]; 
[solObs, qObs, QObs, res] = solvesos(consList, [], options, obsCoeffList); 

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
[sol_free, q_free, Q_free, res] = solvesos(constList_K, (lambda-lambda_free)^2, options, coefList_free); 
% Since K_integral and K_differential share coefficients list, below clean
% statements give the integrated controller. integrableK is the
% differential version, and integralK is the integrated version. 
integrableK = clean(replace(K_differential, coefList_free, double(coefList_free)), 1e-7); 
integralK = clean(replace(K_integral, coefList_free, double(coefList_free)), 1e-7); 
sdisplay(integrableK); 
sdisplay(integralK); 

% %% Second method
% 
% [kappa, kappaCoeff, ~] = polynomial([x1, x2], max(degree(K, [x1,x2])) + 3, 1); 
% dummy = sdpvar(2,1); 
% Qs = sdpvar(2,2); 
% kappaDifferential = jacobian(kappa, [x1 x2]); 
% sosProb = sos( dummy'*(B*kappaDifferential*Wd + Wd*kappaDifferential'*B' - verifiedrho*(B*B') - Qs)*dummy ); 
% % sosProb2 = sos( slack2 - dummy'*(B*kappaDifferential*Wd + Wd*kappaDifferential'*B' - verifiedrho*(B*B') )*dummy ); 
% [solKappa, qKappa, QKappa, res] = solvesos([sosProb; Qs >= 0.01*eye(2); Qs <= 0.5*eye(2) ], [], options, [kappaCoeff; Qs(:)]); 
% 
% kappa = clean(replace(kappa, kappaCoeff, double(kappaCoeff)), 1e-7); 
% % For convenience, I copy and pasted result of sdisplay(integralK) in
% % MGmodel. The last argument is ustar. Because of nonlinearity, for every
% % ustar there is associated equilibrium. 

%% Observer design 
disp('Solving Integrable Observer'); 

Md = double(M); 
[observer, observerCoeff, ~] = polynomial( [x1, x2], max(degree(rhoObs, [x1, x2])) + 2, 1); 
L = jacobian(observer, [x1, x2])'; 
lambdaFreeObs = sdpvar(1); 
M_L_constr = sos(-delta'*(Md*F + F'*Md - Md*L*C - C'*L'*Md + 2*lambdaFreeObs*Md)*delta); 
coeffListObsInt = [observerCoeff(:); lambdaFreeObs];
[solIntObs, qIntObs, QIntObs, ~] = solvesos(M_L_constr, -lambdaFreeObs, options, coeffListObsInt); 
observer = clean(replace(observer, observerCoeff, double(observerCoeff) ), 1e-7); 
sdisplay(observer); 

MGmodelToZero = @(t_anon, x_anon)(MGmodel(t_anon, x_anon, 0)); 
MGmodelToRef = @(t_anon, x_anon)(MGmodel(t_anon, x_anon, 1)); 
[t_to_zero, x_to_zero] = ode45(MGmodelToZero, [0 5], [2 0]');
[t_to_ref, x_to_ref] = ode45(MGmodelToRef, [0 5], [0 0]'); 
figure(1); 
plot(t_to_zero, x_to_zero(:,1)); 
title('Stabilizing Controller - State 1'); 
grid on; 
figure(2); 
plot(t_to_zero, x_to_zero(:, 2)); 
title('Stabilizing Controller - State 2');  
grid on; 
figure(3);
plot(t_to_ref, x_to_ref(:,1)); 
title('Reference Controller to [x1, x2] = [20, 20] - State 1');  
grid on; 
figure(4);
plot(t_to_ref, x_to_ref(:,2)); 
title('Reference Controller to  [x1, x2] = [20, 20] - State 2');  
grid on; 