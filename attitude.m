%% 

clc;
clear; 
close all; 

%% Add MOSEK and YALMIP Path : Please modify appropriately. 

addpath '/home/lkm1321/mosek/7/toolbox/r2013a'; 
addpath(genpath('/home/lkm1321/yalmip')); 

Isq = [2, 0, 0;
           0.1, 3, 0;
           -0.2, 0.3, 5]; 

 I = Isq *Isq';  

 Iwheel = diag([0.1, 0.1, 0.1]); 
 
quat = sdpvar(4,1); 
omega = sdpvar(3,1); 
omegaWheel = sdpvar(3,1); 

x = [quat; omega; omegaWheel]; 

f = [0.5 * quatmultiply([omega; 0]',quat')';
      I\( cross(I*omega, omega) + cross(Iwheel*omegaWheel, omega)) ;
      zeros(3,1)]; 
      
B =  [zeros(4, 3);
         -I\(Iwheel); 
         eye(3)]; 
     
A = jacobian(f, [quat; omega; omegaWheel]); 

% Defining variables. Note YalmipPolynomial class wrapper. 

rho = YalmipPolynomial([x], 2, 0); 
rho2 = YalmipPolynomial([x], 2, 0); 

% Defining matrix variables. Note YalmipPolynomialMatrix. 

W = sdpvar(length(x), length(x)); 
Wdot = zeros(length(x), length(x)); 

% W = YalmipPolynomialMatrix(x, 2, 0, length(x), length(x)); 
% Wdot = W.lieDerivative(f); 

%W = getPolynomialMatrix(x, 2, 0, length(x), length(x)); 
% Wdot = getMatrixDot(W, f, x); 

% w11 = YalmipPolynomial([x], 2, 0);
% w12 = YalmipPolynomial([x], 2, 0); 
% w22 = YalmipPolynomial([x], 2, 0); 
% W = [w11.polyRep, w12.polyRep; 
%         w12.polyRep, w22.polyRep]; 
% 
%  w11dot = jacobian(W(1,1), [x1, x2]) * f;
%  w12dot = jacobian(W(1,2), [x1, x2]) * f;
%  w22dot = jacobian(W(2,2), [x1, x2]) * f;
% 
%  Wdot = [w11dot, w12dot; 
%                w12dot, w22dot];

lambda = 0.25;
ep = 1e-6;                                  % Defining tolerences
upper = 10; 
lower = 0.1;                                % Upper and lower bounds on W. 

delta = sdpvar(length(x),1);                        % Dummy variable
           
LMI = -Wdot + A*W + W*A' - rho.polyRep * (B*B') - rho2.polyRep *(f*f')+2*lambda*W; 

constraints = [sos(-delta'*LMI*delta);
                        sos(delta' *(W  -  lower*eye(size(W)))*delta);
                        sos(delta' *(upper*eye(size(W)) - W)*delta);
                        sos(rho.polyRep)
                        sos(rho2.polyRep)];
                    
coeffList = [W
                    rho.coeffList
                    rho2.coeffList];
                
[sol, q, Q, res] = solvesos(constraints, [], options, coeffList);
% W = clean(replace(W, [w11.coeffList; w12.coeffList; w22.coeffList], double([w11.coeffList; w12.coeffList; w22.coeffList]))); 
% Fix W and rho. 
W = double(W); 
rho.polyRep = clean(replace(rho.polyRep, rho.coeffList, double(rho.coeffList)), 1E-7);  
rho2.polyRep = clean(replace(rho2.polyRep, rho2.coeffList, double(rho2.coeffList)), 1E-7); 
%% Finding integrable controller. 
controller = YalmipPolynomial(x, 2, 1); 
% virtualController = YalmipPolynomial([x1, x2], 2, 1); 
differential = jacobian(controller.polyRep, x);
% virtualDifferential = jacobian(virtualController.polyRep, [x1, x2]); 
K = [differential]; 

freeLambda = sdpvar(1); 

LMI2 = (-Wdot + A*W + W*A' - B*K*W - W*K'*B' - rho2.polyRep*(f*f') + 2*freeLambda*W );

intConstraints = [sos(-delta'*LMI2*delta);  
                              (freeLambda >= 0)];

intCoeffList = [controller.coeffList; 
                        freeLambda]; 

[sol_free, q_free, Q_free, res] = solvesos(intConstraints, -freeLambda, options, intCoeffList);

verifiedController = clean(replace(controller.polyRep, controller.coeffList, double(controller.coeffList)), 1E-7); 
verifiedDifferential = clean(replace(differential, controller.coeffList, double(controller.coeffList)), 1E-7); 

     