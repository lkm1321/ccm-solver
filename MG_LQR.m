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

A = clean(replace(F, [x1 x2], [0 0], 1e-7));

Q = 4*eye(2);
R = 4; 

K = lqr(A, B, Q, R); 
