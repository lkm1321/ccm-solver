%% Van Der Pol Equations 
%-------------------------------------------------------------------------%
% This code finds the contraction metrics for the Van Der Pol system given
% the Lienard substitution 2D system. The equations are taken from
% Guckenheimer paper. Three different oscillators are considered and a
% contraction metric that satisfies all is considered
% Note: the degree of rho for the observer is 4 while the degree for the
% controller is 2.
% Author: Karen Leung 310241847
% -- Last updated 10/6/2014 -- %
%-------------------------------------------------------------------------%
function [W,c1,WO,c2,tout1,xout1,tout2,xout2,t_lead,x_lead,t_back1,x_back1,t_back2,x_back2] = Feedback_VdPSingleControl()
addpath 'c:\Program Files\mosek\7\toolbox\r2013a'

x1 = sdpvar(1);                         % state variable 1
x2 = sdpvar(1);                         % state variable 2
ep = 1e-6;                              % tolerence
delta = sdpvar(2,1);                    % quadratic variable
B = [0;1];                              % B controller matrix
C = [1 0];                              % C observer matrix
x = [x1 x2]';                           % state space
mu = 5;                                 % Van Der Pol parameter
dmu1 = 1.4;
dmu2 = 1.8;
xdot = mu * (x2 - (x1^3/3 - x1));       % state equation 1    
ydot = - x1 / mu;                       % state equation 2
f1 = [xdot ; ydot];                     % state equation     
F1 = jacobian(f1,[x1 x2]);              % jacobian matrix

xdot = mu*dmu2 * (x2 - (x1^3/3 - x1));  % state equation 1    
ydot = - x1 / (mu*dmu2);                % state equation 2
f2 = [xdot ; ydot];                     % state equation     
F2 = jacobian(f2,[x1 x2]);              % jacobian matrix


%% Control Contraction Metric set up
r_rho = 2;                              % degree of polynomial for rho
lambda = 12.23;                         % Convergence rate
W = sdpvar(2,2);                        % polynomials for W matrix

[rho , c1 , v1] = polynomial([x1 x2],r_rho,0);  % rho polynomial
H1 = (F1*W+W*F1'-rho*B*B'+2*lambda*W);  % controller contraction matrix
Wdotconstr1 = sos(-delta'*H1*delta);    % inequality SOS
H2 = (F2*W+W*F2'-rho*B*B'+2*lambda*W);  % controller contraction matrix
Wdotconstr2 = sos(-delta'*H2*delta);    % inequality SOS
lower = 0.1;                            % lower bound
Wconstr = [W>=lower*eye(2)];
upper = 5;                              % upper bound
Wconstr2 = [W<=upper*eye(2)];
rhoconstr = sos(rho);                   % rho SOS

% Solving the control contraction metric
options = sdpsettings('solver','mosek');
consList = [Wdotconstr1,Wdotconstr2, Wconstr, Wconstr2,rhoconstr];
coefList=[W(:); c1];
% search for SOS certificate
[sol, q, Q, res] = solvesos(consList, [], options, coefList);

verifiedrho = clean(replace(rho, coefList, double(coefList)), 1e-7);
sdisplay(verifiedrho)

Wd = double(W)

%% Observer Contraction set up
r_rho = 4;                              % degree of polynomial for rho
% lambda = 45.5;                        % convergence rate
lambda = 35.6;
WO = sdpvar(2,2);                       % polynomials for W matrix

[rho , c2 , v2] = polynomial([x1 x2],r_rho,0);  % rho polynomial
H1=(F1'*WO+WO*F1-rho*C'*C+2*lambda*WO); % Observer Contraction Matrix
Wdotconstr1 = sos(-delta'*H1*delta);    % WO is SOS
H2=(F2'*WO+WO*F2-rho*C'*C+2*lambda*WO); % Observer Contraction Matrix
Wdotconstr2 = sos(-delta'*H2*delta);    % WO is SOS
lower = 0.01;                           % lower bound
Wconstr = [WO>=lower*eye(2)];
rhoconstr = sos(rho);                   % rho is SOS

% Solving for observer contraction metric
options = sdpsettings('solver','mosek');
consList = [Wdotconstr1,Wdotconstr2, Wconstr, rhoconstr];
coefList=[c2; WO(:)];
% search for SOS certificate
[sol, q, Q, res] = solvesos(consList, [], options, coefList);

verifiedrho = clean(replace(rho, coefList, double(coefList)), 1e-7);
sdisplay(verifiedrho)

WdO = double(WO)



%% Constructing controller and observer and integrating
x0 = [0.7;0.6;0.5;0.5];                 % current and observed state 
ustar = 0;                              % initial control

rholistC = double(c1);                  % coefficients of controller rho
rholistO = double(c2);                  % coefficients of observer rho 

t_max = 20;                             % Max time of simulation
lead = @(t,x) [mu*dmu1*(x(2)-(x(1)^3/3-x(1)));-x(1)/(mu*dmu1)]; % Desired trajectory 
[t_lead,x_lead] = ode45(lead,[0,t_max],x0(1:2)*2*dmu1);
t_on = 10;                              % Time when controller is turned on

back1 = @(t,x) [mu*(x(2)-(x(1)^3/3-x(1)));-x(1)/(mu)]; % Before controller is turned on
[t_back1,x_back1] = ode45(back1,[0,t_on],x0(1:2));
x_on1 = x_back1(end,:);

back2 = @(t,x) [mu*dmu2*(x(2)-(x(1)^3/3-x(1)));-x(1)/(mu*dmu2)]; 
[t_back2,x_back2] = ode45(back2,[0,t_on],x0(1:2)*2*dmu2);
x_on2 = x_back2(end,:);

tic
[tout1,xout1]= ode45(@(t,xcurr) Construct_VdPSingleControl(t,xcurr,Wd,WdO,rholistC,rholistO, ustar,mu,t_lead,x_lead), [t_on t_max], [x_on1';x_on1'*1.5]);
[tout2,xout2]= ode45(@(t,xcurr) Construct_VdPSingleControl(t,xcurr,Wd,WdO,rholistC,rholistO, ustar,mu,t_lead,x_lead), [t_on t_max], [x_on2';x_on2'*1.5]);
toc

