%% xdot = Ax 
%% y = Cx

%A = [-1, 1;
%      0, 1]; 
clear;
clc;
close all;


C = [2, 1]; 
Cperp = [-1, 2]'; 

V = [Cperp, C']; 
D = diag([1,-2]);

A = V*D* (V^-1);

numSamples = 20; 
numTSamples = 100; 
tEnd = 1; 

tRange = linspace(0, tEnd, numTSamples); 

x1Range = linspace(-10, 10, numSamples); 
x2Range = linspace(-10, 10, numSamples); 

 

x1Dot = zeros(length(x1Range), length(x2Range) );
x2Dot = zeros(length(x1Range), length(x2Range) ); 

x1Response = zeros( length(tRange), length(x1Range), length(x2Range) ); 
x2Response = x1Response; 
obsResponse = x1Response; 

xdot = @(t,x)(A*x); 
obs = @(x)(C*x);

t = 0; 

for i = 1:length(x1Dot)
    for j = 1:length(x2Dot)
        x = [x1Range(i);x2Range(j)]; 
        xDot = xdot(t,x); 
        x1Dot(i,j) = xDot(1);
        x2Dot(i,j) = xDot(2); 
        
        for timeIdx = 1:length(tRange)
            xResponse = expm(A*tRange(timeIdx))*x; 
            obsResponse(timeIdx, i, j) = obs(x); 
            x1Response(timeIdx, i, j) = xResponse(1);
            x2Response(timeIdx, i, j) = xResponse(2);
        end
        
    end
end

[x1, x2] = meshgrid(x1Range, x2Range);
 
figure(1)
quiver(x1,x2,x1Dot, x2Dot); 

figure(2)
hold on; 
for i = 1:length(x1Range)
    for j = 1:length(x2Range)
        plot(x1Response(:,i,j), x2Response(:,i,j)); 
    end
end
hold off; 

figure(3)
hold on; 
for i = 1:length(x1Range)
    for j = 1:length(x2Range)
        plot(obsResponse(:,i,j)); 
    end
end
hold off; 




