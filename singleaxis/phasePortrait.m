clear;
clc;
close all;
clear controller; 
clear f;

deFunc = @(t,x)(fnew(x,controllerNew(x))); 

numSample = 100; 

thetaRange = linspace(-2*pi, 2*pi, numSample); 
omegaRange = linspace(-pi, pi, numSample); 

thetaDot = zeros(length(thetaRange), length(omegaRange)); 
omegaDot = zeros(length(thetaRange), length(omegaRange)); 

t = 0; 

for i = 1:length(thetaRange)
    for j = 1:length(omegaRange)
        
        x = [thetaRange(i); omegaRange(j)]; 
        xDot = deFunc(t,x); 
        thetaDot(i,j) = xDot(1);
        omegaDot(i,j) = xDot(2); 
        
    end 
end

[x,y] = meshgrid(thetaRange, omegaRange); 

figure(1)
quiver(x,y, thetaDot, omegaDot); 

figure(2)
surf(x,y,thetaDot);
xlabel('Theta'); 
ylabel('Omega'); 

figure(3)
surf(x,y,omegaDot);
xlabel('Theta');
ylabel('Omega'); 



