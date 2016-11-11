clc;
clear;
close all; 

% Uncontrolled LV system. 

uncontrolledInitial1 = [1.5, 0.5]'; 
uncontrolledInitial2 = [2, 2]'; 
uncontrolledInitial3 = [3, 4]'; 


[tu1, xu1] = ode45(@LVuncontrolled, [0, 10], uncontrolledInitial1); 
[tu2, xu2] = ode45(@LVuncontrolled, [0, 10], uncontrolledInitial2); 
[tu3, xu3] = ode45(@LVuncontrolled, [0, 10], uncontrolledInitial3); 

% Controlled LV system. Stabilizes to a point

controlledInitial1 = [1, 1]';  
controlledInitial2 = [1.5, 0.5]';
controlledInitial3 = [2, 2]';  

[tc1, xc1] = ode45(@LVcontrolled, [0, 3.5], controlledInitial1); 
[tc2, xc2] = ode45(@LVcontrolled, [0, 3.5], controlledInitial2);
[tc3, xc3] = ode45(@LVcontrolled, [0, 3.5], controlledInitial3); 

% Reference tracking

ref1 = 2; 
ref2 = 3; 
ref3 = 4; 

refModel1 = @(t, x)(LVreference(t, x, ref1));
refModel2 = @(t, x)(LVreference(t, x, ref2));
refModel3 = @(t, x)(LVreference(t, x, ref3)); 

[tr1, xr1] = ode45(refModel1, [0, 10], [1, 1]'); 
[tr2, xr2] = ode45(refModel2, [0, 10], [1, 1]'); 
[tr3, xr3] = ode45(refModel3, [0, 10], [1, 1]'); 


% Natural Manifold Stabilization. 

naturalFinal1 = [1, 1, 2, 2]'; 
naturalFinal2 = [1, 1, 2, 3]'; 
naturalFinal3 = [1, 1, 3, 2]'; 

[tn1, xn1] = ode45(@LVstabilized, [0, 100], naturalFinal1);
[tn2, xn2] = ode45(@LVstabilized, [0, 100], naturalFinal2); 
[tn3, xn3] = ode45(@LVstabilized, [0, 100], naturalFinal3); 

% Circle Stabilization 

circleFinal1 = [1,1, 1.25, 1]'; 
circleFinal2 = [1,1, 1.5, 1]'; 
circleFinal3 = [1,1, 2, 1]'; 

[ti1, xi1] = ode45(@LVsynchronizeToCircle, [0, 100], circleFinal1); 
[ti2, xi2] = ode45(@LVsynchronizeToCircle, [0, 100], circleFinal2);
[ti3, xi3] = ode45(@LVsynchronizeToCircle, [0, 100], circleFinal3);

% Circle stabilization from different initial points

circ2F1 = [0.5, 3, 1.5, 1]'; 
circ2F2 = [2, 3, 1.5, 1]'; 
circ2F3 = [2, 2, 1.5, 1]'; 

[tF1, xF1] = ode45(@LVsynchronizeToCircle, [0, 200], circ2F1); 
[tF2, xF2] = ode45(@LVsynchronizeToCircle, [0, 200], circ2F2);
[tF3, xF3] = ode45(@LVsynchronizeToCircle, [0, 200], circ2F3);

% Plots
figure; 
hold on;
plot(xu1(:,1), xu1(:,2));
plot(xu2(:,1), xu2(:,2));
plot(xu3(:,1), xu3(:,2)); 
hold off; 
title('Uncontrolled Limit Cycles'); 
xlabel('x1'); 
ylabel('x2'); 
grid on; 

figure;
hold on; 
plot(xc1(:,1), xc1(:,2)); 
plot(xc2(:,1), xc2(:,2)); 
plot(xc3(:,1), xc3(:,2)); 
hold off; 
title('Stabilization to zero'); 
xlabel('x1');
ylabel('x2'); 
grid on; 

figure;
subplot(2, 1, 1); 
hold on; 
plot(tc1(:,1), xc1(:,1)); 
plot(tc2(:,1), xc2(:,1)); 
plot(tc3(:,1), xc3(:,1)); 
hold off; 
title('Stabilization to zero'); 
% xlabel('t');
ylabel('x1'); 
grid on; 
subplot(2,1, 2); 
hold on; 
plot(tc1(:,1), xc1(:,2)); 
plot(tc2(:,1), xc2(:,2)); 
plot(tc3(:,1), xc3(:,2)); 
hold off; 
% title('Stabilization to natural manifold'); 
xlabel('t');
ylabel('x2'); 
grid on;

figure;
subplot(2, 1, 1); 
hold on; 
plot(tr1(:,1), xr1(:,1)); 
plot(tr2(:,1), xr2(:,1)); 
plot(tr3(:,1), xr3(:,1)); 
hold off; 
title('Stabilization to reference point'); 
% xlabel('t');
ylabel('x1'); 
grid on; 
subplot(2,1, 2); 
hold on; 
plot(tr1(:,1), xr1(:,2)); 
plot(tr2(:,1), xr2(:,2)); 
plot(tr3(:,1), xr3(:,2)); 
hold off; 
% title('Stabilization to natural manifold'); 
xlabel('t');
ylabel('x2'); 
grid on;


figure; 
hold on; 
plot(xn1(:,1), xn1(:,2)); 
plot(xn2(:,1), xn2(:,2)); 
plot(xn3(:,1), xn3(:,2)); 
hold off; 
title('Stabilization to natural manifold'); 
xlabel('x1');
ylabel('x2'); 
grid on; 

figure; 
hold on; 
plot(xi1(:,1), xi1(:,2)); 
plot(xi2(:,1), xi2(:,2)); 
plot(xi3(:,1), xi3(:,2)); 
hold off; 
title('Stabilization to Different Circles from Equilibrium'); 
xlabel('x1');
ylabel('x2'); 
grid on; 
axis equal; 

figure; 
hold on;
plot(xF1(:,1), xF1(:,2));
plot(xF2(:,1), xF2(:,2));
plot(xF3(:,1), xF3(:,2)); 
hold off; 
title('Stabilization to Circle from Different Initial Conditions'); 
xlabel('x1'); 
ylabel('x2'); 
grid on; 

