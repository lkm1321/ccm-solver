%% Contraction Analysis via SOS programming Moore Greitzer Model
%-------------------------------------------------------------------------%
% This function generates the algorithm presented by Aylward, Parrilo and
% Slotine for testing if a system is contracting. Here, the Moore-Greitzer
% model is tested. The metric is modelled as a 4th order polynomial matrix
% and through a bisection process, the bounds for the rate of convergence
% can be found. 
% Author: Karen Leung 310241847
% -- Last updated 9/25/2014 -- %
%-------------------------------------------------------------------------%
function [MM] = Contraction_MG()
addpath 'c:\Program Files\mosek\7\toolbox\r2013a'
% phi = x
% psi = y
%% Setting up variables
sdpvar x y
b = 0.78;                         % Convergence rate
Y = sdpvar(2,1);                  % Dummy variable  
d = 4;                            % Degree of metric  

%% The system equation and its jacobian
xdot = -y - 3/2*x^2 - 0.5*x^3;
ydot =  3*x - y;
f = [xdot ; ydot];                % Setting up the EOM
dfdx = jacobian(f , [x , y]);     % Calculating Jacobian

%% Setting up M and Mdot matrix
variables = monolist([x,y],d);      
variablesdot = variables;
[n , ~] = size(variables);

A = sdpvar(n,1);
B = sdpvar(n,1);
D = sdpvar(n,1);

M = sdpvar(2);
Mdot = sdpvar(2);
M(1,1) = A'*variables;
M(1,2) = B'*variables;
M(2,1) = B'*variables;
M(2,2) = D'*variables;

for i = 1:n
    variablesdot(i) = jacobian(variables(i),[x , y])*[xdot ; ydot];
end

Mdot(1,1) = A'*variablesdot;
Mdot(1,2) = B'*variablesdot;
Mdot(2,1) = B'*variablesdot;
Mdot(2,2) = D'*variablesdot;
    
%% Setting up constraint
R = dfdx'*M + M*dfdx + Mdot + b*M;              % LMI constraint
con1 = sos(Y'*M*Y);
con2 = sos(-Y'*R*Y);
constraints = [con1 , con2];
coefficients = [A ; B ; D ];
options = sdpsettings('solver','mosek');
F = solvesos(constraints, [], options, coefficients);
MM = clean(replace(M, coefficients, double(coefficients)), 1e-4);
fprintf('The contraction metric is:\n')
sdisplay(MM(1,1))
sdisplay(MM(1,2))
sdisplay(MM(2,1))
sdisplay(MM(2,2))





