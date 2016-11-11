clear;
clc;
close all;

xstar1 = 0; 
xstar2 = 2; 
xstar3 = 4; 
xstar4 = 6;  

MGUncontrolled = @(t, x) (  [-0.5*x(1)^3-1.5*x(1)^2-x(2); x(1)] ); 
MGStabilize = @(t, x)(MGmodel(t, x, xstar1)); 
MGTrack1 = @(t, x)(MGmodel(t, x, xstar2)); 
MGTrack2 = @(t, x)(MGmodel(t, x, xstar3)); 
MGTrack3 = @(t, x)(MGmodel(t, x, xstar4)); 

[t0, x0] = ode45(MGUncontrolled, [0, 100], [10 5]'); 
[tPid, xPid] = ode45(@(t,x)MGpd(t,x,0), [0, 10], [2 2]'); 
% Stabilizing
[t1, x1] = ode45(MGStabilize, [0, 1.5], [10 5]'); 
[t2, x2] = ode45(MGStabilize, [0, 1.5], [5 10]'); 
[t3, x3] = ode45(MGStabilize, [0, 1.5], [5 5]'); 
% Tracking
[t4, x4] = ode45(MGTrack1, [0, 1.5], [0 0]');
[t5, x5] = ode45(MGTrack2, [0, 1.5], [0 0]');
[t6, x6] = ode45(MGTrack3, [0, 1.5], [0 0]'); 
% Observing
% [tObs1, xObs1] = ode45(@(t,x)(MGobserverModel(t,x,2)), [0, 2.5], [0 0 10 10]'); 
% Output Feedback
% [t, x] = ode45(@(t,x)(MGoutputFeedbackSecond(t,x,0)), [0, 2], [0 0 10 10]');
% I somehow messed up the order when I was writing the thesis. Better to
% change code. 

figure; 
hold on; 
plot(t1, x1(:,1));
plot(t1, x1(:,2), 'r');
hold off; 
title('Stabilizing controller');
grid on;
xlabel('Time (s)');
legend('x1', 'x2'); 

figure;
hold on;
plot(t1, x1(:,1)); 
plot(t2, x2(:,1));
plot(t3, x3(:,1)); 
hold off;
title('Stabilizing Controller - State 2'); 
grid on; 
xlabel('Time (s)');
ylabel('x2'); 

figure;
hold on;
plot(t1, x1(:,2)); 
plot(t2, x2(:,2));
plot(t3, x3(:,2)); 
hold off;
title('Stabilizing Controller - State 1'); 
grid on; 
xlabel('Time (s)');
ylabel('x1'); 

figure;
hold on;
plot(tPid, xPid(:,1)); 
plot(tPid, xPid(:,2)); 
hold off;
title('Why not PID?'); 
grid on;
xlabel('Time (s)'); 
legend('x1', 'x2');

figure;
hold on; 
plot(t0, x0(:,1)); 
plot(t0, x0(:,2)); 
hold off;
legend('State x2', 'State x1'); 
title('Uncontrolled MG Dynamics');
grid on;
xlabel('Time (s)'); 

figure; 
hold on;
plot(t4, x4(:, 1));
plot(t5, x5(:, 1));
plot(t6, x6(:, 1)); 
hold off; 
title('Tracking Controller - x2');
grid on;
xlabel('Time (s)'); 

figure; 
hold on;
plot(t4, x4(:, 2));
plot(t5, x5(:, 2));
plot(t6, x6(:, 2)); 
hold off; 
title('Tracking Controller - x1');
grid on;
xlabel('Time (s)'); 

figure; 
hold on; 
plot(tObs1, xObs1(:,1), tObs1, xObs1(:,3));
plot(tObs1, xObs1(:,2), tObs1, xObs1(:,4)); 
hold off; 
title('State Estimator'); 
legend('xhat1', 'xtrue1', 'xhat2', 'xtrue2'); 
grid on; 
xlabel('Time (s)'); 

rms = sqrt( sum(0.2*randn(size(xObs1(:,1),1),1).^2)/size(xObs1,1) )*ones(size(xObs1,1),1); 

figure; 
hold on; 
plot(tObs1(50:end,1), sqrt((xObs1(50:end,4) - xObs1(50:end, 2)).^2), tObs1(50:end,1), rms(50:end)); 
hold off;
title('Noise Rejection'); 
legend('Absolute Error in x2', 'RMS of noise'); 
axis([0.5, max(tObs1), 0, max(rms)+1]); 
grid on ; 
xlabel('Time(s)');

figure; 
hold on;
plot(t, x(:,1));
plot(t, x(:, 2));
plot(t, x(:,3));
plot(t, x(:, 4));
hold off;
title('Output Feedback'); 
legend('xhat1', 'xhat2', 'x1', 'x2'); 
grid on;
xlabel('Time (s)');
