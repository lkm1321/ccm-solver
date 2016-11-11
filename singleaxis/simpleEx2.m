C = [2, 1];
Cperp = [-1, 2]; 

D = diag([1, -1]); 
V = [C', Cperp']; 
A = V*D*V^-1; 

sys = ss(A, [], C, []); 
sys2 = ss(A, [], Cperp, []); 

x0 = [1; 1]; 

figure, initial(sys, x0, 10); 

figure, initial(sys2, x0, 10); 