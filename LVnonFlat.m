clear;
clc;
close all; 

%% Add MOSEK and YALMIP Path : Please modify appropriately. 

addpath '/home/lkm1321/mosek/7/toolbox/r2013a'; 
addpath(genpath('/home/lkm1321/yalmip')); 

%% Finding CCM. Parts taken from the works of Karen Leung and Max Revay. 

x1 = sdpvar(1);                             % Defining variables
x2 = sdpvar(1);                             % Defining variables

x = [x1 x2]';                               % State space


f = [x1-x1*x2;
      -x2+x1*x2];             % Dynamical system
B = [1;0];                                  % Control matrix
C = [0 1];                                  % Observation matrix

A = jacobian(f,[x1 x2]);                    % Jacobian matrix

% Defining variables. Note YalmipPolynomial class wrapper. 

rho = YalmipPolynomial([x1, x2], 2, 0); 
rho2 = YalmipPolynomial([x1, x2], 2, 0); 
w11 = YalmipPolynomial([x1, x2], 2, 0);
w12 = YalmipPolynomial([x1, x2], 2, 0); 
w22 = YalmipPolynomial([x1, x2], 2, 0); 
W = [w11.polyRep, w12.polyRep; 
        w12.polyRep, w22.polyRep]; 

 w11dot = jacobian(W(1,1), [x1, x2]) * f;
 w12dot = jacobian(W(1,2), [x1, x2]) * f;
 w22dot = jacobian(W(2,2), [x1, x2]) * f;

 Wdot = [w11dot, w12dot; 
               w12dot, w22dot];

lambda = 0.25;
ep = 1e-6;                                  % Defining tolerences
upper = 10; 
lower = 0.1;                                % Upper and lower bounds on W. 

delta = sdpvar(2,1);                        % Dummy variable
           
LMI = -Wdot + A*W + W*A' - rho.polyRep * (B*B') - rho2.polyRep *(f*f')+2*lambda*W; 

constraints = [sos(-delta'*LMI*delta);
                        sos(delta' *(W  -  lower*eye(size(W)))*delta);
                        sos(delta' *(upper*eye(size(W)) - W)*delta);
                        sos(rho.polyRep)
                        sos(rho2.polyRep)];
                    
coeffList = [w11.coeffList; 
                    w12.coeffList;
                    w22.coeffList;
                    rho.coeffList
                    rho2.coeffList]; 

options = sdpsettings('solver','mosek');

[sol, q, Q, res] = solvesos(constraints, [], options, coeffList);
W = clean(replace(W, [w11.coeffList; w12.coeffList; w22.coeffList], double([w11.coeffList; w12.coeffList; w22.coeffList]))); 
rho.polyRep = clean(replace(rho.polyRep, rho.coeffList, double(rho.coeffList)), 1E-7);  
rho2.polyRep = clean(replace(rho2.polyRep, rho2.coeffList, double(rho2.coeffList)), 1E-7); 
M = [W(2,2), -W(1,2);
        -W(1,2),  W(1,1)]/(W(1,1)*W(2,2) - W(1,2)^2); 
 K= rho.polyRep*B'*M; 

%% Finding integrable controller. 
controller = YalmipPolynomial([x1, x2], 3, 1); 
differential = jacobian(controller.polyRep, [x1 x2]); 
freeLambda = sdpvar(1); 

LMI2 = (A*W + W*A' - B*differential*W - W*differential'*B' - rho2.polyRep*(f*f') + 2*freeLambda*W );

intConstraints = [sos(-delta'*LMI2*delta);  
                              freeLambda >= 0];

intCoeffList = [controller.coeffList;
                        freeLambda]; 

[sol_free, q_free, Q_free, res] = solvesos(intConstraints, (lambda-freeLambda)^2, options, intCoeffList);

verifiedController = clean(replace(controller.polyRep, controller.coeffList, double(controller.coeffList)), 1E-7); 
verifiedDifferential = clean(replace(differential, controller.coeffList, double(controller.coeffList)), 1E-7); 


                          
                              





    
    

