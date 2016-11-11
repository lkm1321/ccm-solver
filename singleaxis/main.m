clear;
clc;
close all;
clear controller; 
clear f; 

dt = 0.1; 
t_end = 2000; 
t = 0; 
i = 1; 
x = zeros(2,ceil(t_end/dt)+1); 
u = zeros(1,ceil(t_end/dt));
x(1,1) = 100*pi/180; 
global k;
if (isempty(k)) 
   k = 500*10^-6/0.005;  
end

while t<t_end
    u(:,i) = controller(x(1,i),x(2,i),k*3);
    x(:,i+1) = RK4(@f,x(:,i),dt); 
    t = t +dt;
    i = i+1; 
end
c = 1;
figure(1);
plot(linspace(0,t_end,i),x(1,:)*180/pi);
title('Angle');
xlabel('Epoch (s)');
ylabel('Yaw (deg)'); 
figure(2);
plot(linspace(0,t_end,i),x(2,:));
title('Velocity');
xlabel('Epoch (s)'); 
ylabel('Angular Velocity (rad/s)'); 
figure(3);
plot(linspace(0,t_end,i-1),u);
title('Control');
xlabel('Epoch (s)'); 
ylabel('Dipole moment (Am^2)'); 
