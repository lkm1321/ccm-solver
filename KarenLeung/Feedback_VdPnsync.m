%% Van Der Pol Equations 
function [W,WO,c1,c2,tout,xout,t_back,x_back,t_lead,x_lead] = Feedback_VdPnsync(dmu)
%-------------------------------------------------------------------------%
% This code synchronises n (typically 3) systems of the unforced Van Der
% Pol oscillator. The convergence rate used is for the mu = 5 case.
% Note: the degree of rho for the observer is 4 while the degree for the
% controller is 2.
% Input: dmu - a vector with 1 in the first entry and a scaling factor for
% all the other oscillators of interest
% Author: Karen Leung 310241847
% -- Last updated 10/5/2014 -- %
%-------------------------------------------------------------------------%
addpath 'c:\Program Files\mosek\7\toolbox\r2013a'
close all
clc

sys = length(dmu);

x1 = sdpvar(1);                         % state variable 1
x2 = sdpvar(1);                         % state variable 2
ep = 1e-6;                              % tolerence
delta = sdpvar(2,1);                    % quadratic variable

x = [x1 x2]';                           % state space
mu = 5*dmu;

B = [0;1];                              % B controller matrix
C = [1 0];                              % C observer matrix

f = cell(sys,1);
F = cell(sys,1);

for i = 1:sys
    xdot = mu(i) * (x2 - (x1^3/3 - x1));       % state equation 1    
    ydot = - x1 / mu(i);                       % state equation 2
    f{i} = [xdot ; ydot];                       % state equation     
    F{i} = jacobian(f{i},[x1,x2]);
end

%% Control Contraction Metric set up

W = cell(sys,1);
c1 = cell(sys,1);
WO = cell(sys,1);
c2 = cell(sys,1);

for i = 1:sys
    r_rho = 2;                              % degree of polynomial for rho
    % lambda = 12.18;                        % convergence rate
    lambda = 12.18;
    W{i} = sdpvar(2,2);                        % polynomials for W matrix
    [rho , c1{i} , v1] = polynomial([x1 x2],r_rho,0);  % rho polynomial
    H = (F{i}*W{i}+W{i}*F{i}'-rho*B*B'+2*lambda*W{i});     % controller contraction matrix
    Wdotconstr = sos(-delta'*H*delta);      % inequality SOS
    lower = 0.1;                            % lower bound
    Wconstr = [W{i}>=lower*eye(2)];
    upper = 5;                              % upper bound
    Wconstr2 = [W{i}<=upper*eye(2)];
    rhoconstr = sos(rho);                   % rho SOS

    % Solving the control contraction metric
    options = sdpsettings('solver','mosek');
    consList = [Wdotconstr, Wconstr, Wconstr2,rhoconstr];
    coefList=[W{i}(:); c1{i}];
    % search for SOS certificate
    [sol, q, Q, res] = solvesos(consList, [], options, coefList);

    verifiedrho = clean(replace(rho, coefList, double(coefList)), 1e-7);
    sdisplay(verifiedrho)
    c1{i} = double(c1{i});
    W{i} = double(W{i})
% Observer Contraction set up
    r_rho = 4;                              % degree of polynomial for rho
    % lambda = 45.5;                        % convergence rate
    lambda = 45.5;
    WO{i} = sdpvar(2,2);                       % polynomials for W matrix
    [rho , c2{i} , v2] = polynomial([x1 x2],r_rho,0);  % rho polynomial
    H=(F{i}'*WO{i}+WO{i}*F{i}-rho*C'*C+2*lambda*WO{i});    % Observer Contraction Matrix
    Wdotconstr = sos(-delta'*H*delta);      % WO is SOS
    lower = 0.01;                           % lower bound
    Wconstr = [WO{i}>=lower*eye(2)];
    rhoconstr = sos(rho);                   % rho is SOS

    % Solving for observer contraction metric
    options = sdpsettings('solver','mosek');
    consList = [Wdotconstr, Wconstr, rhoconstr];
    coefList=[c2{i}; WO{i}(:)];
    % search for SOS certificate
    [sol, q, Q, res] = solvesos(consList, [], options, coefList);

    verifiedrho = clean(replace(rho, coefList, double(coefList)), 1e-7);
    sdisplay(verifiedrho)
    c2{i} = double(c2{i});
    WO{i} = double(WO{i})
end


%% Constructing controller and observer and integrating
x0 = [0.7;0.6;0.5;0.5];                 % current and observed state 
ustar = 0;                              % initial control
t_max = 50;
lead = @(t,x) [mu(1)*(x(2)-(x(1)^3/3-x(1)));-x(1)/(mu(1))];  
[t_lead,x_lead] = ode45(lead,[0,t_max],x0(1:2));
t_on = 20;

back = cell(sys,1);
t_back = cell(sys,1);
x_back = cell(sys,1);
x_on = cell(sys,1);
t_out = cell(sys,1);
x_out = cell(sys,1);

for i = 2:sys
    back{i} = @(t,x) [mu(i)*(x(2)-(x(1)^3/3-x(1)));-x(1)/(mu(i))]; 
    [t_back{i},x_back{i}] = ode45(back{i},[0,t_on],x0(1:2)*2*dmu(i));
    x_on{i} = x_back{i}(end,:);
    tic
    [t_out{i},x_out{i}]= ode45(@(t,xcurr) Construct_VanDerPol(t,xcurr,W{i},WO{i},c1{i},c2{i},ustar,mu(i),t_lead,x_lead), [t_on t_max], [x_on{i}';x_on{i}'*1.5]);
    toc    
end


