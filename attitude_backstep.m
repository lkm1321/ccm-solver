%% 

clc;
clear; 
close all; 

%% Add MOSEK and YALMIP Path : Please modify appropriately. 

addpath '/home/lkm1321/mosek/7/toolbox/r2013a'; 
addpath(genpath('/home/lkm1321/yalmip')); 

%% Kinematics

quat = sdpvar(4,1);
omega = sdpvar(3,1); 
x = [quat; omega]; 

f = [0.5 * quatmultiply([omega; 0]',quat')']; 

A = jacobian(f, quat); 
B = jacobian(f, omega); 

G = [-quat(2), -quat(3), -quat(4);
        quat(1), -quat(4), -quat(3); 
        quat(4), quat(1), quat(2); 
        quat(3), quat(2), quat(1)]; 

W = sdpvar(4, 4); 
Wdot = zeros(4, 4); 

rho = YalmipPolynomial( [quat; omega], 2, 0); 
rho2 = YalmipPolynomial( [quat; omega], 2, 0); 

lambda = 0.0;
ep = 1e-6;                                  % Defining tolerences
upper = 100; 
lower = 0.1;                                % Upper and lower bounds on W. 

delta = sdpvar(4,1);                        % Dummy variable
delta2 = sdpvar(3, 1); 
gamma = sdpvar(1); 

LMI = -Wdot + A*W + W*A' -  (B*rho.polyRep*B') - (G*rho2.polyRep*G')+2*lambda*W; 

constraints = [sos(-delta'*LMI*delta);
                        (W >= lower * eye(4)); 
                        (W <= upper * eye(4));
                        sos(rho.polyRep); 
                        sos(rho2.polyRep)];

coeffList = [W(:);
                    rho.coeffList;
                    rho2.coeffList];
                
options = sdpsettings('solver','mosek', 'verbose', 2);
                
[sol, q, Q, res] = solvesos(constraints, [], options, coeffList);                
