clear
clc

x = sdpvar;
y = sdpvar;

%Moore Greitzer model
fx = [-y  - 1.5*x^2 - 0.5*x^3;
       x];
   
 A = jacobian(fx,[x,y]);
 B = [0;1];
 
 monos = monolist([x,y],degree(fx));
 
 %coefficient of our polynomials in W and also decision variables
 w11 = sdpvar(length(monos),1);
 w22 = sdpvar(length(monos),1);
 w12 = sdpvar(length(monos),1);  
 
 % The matrix W is symmetrix with a polynomial in each element.
 W = [w11'*monos, w12'*monos;
      w12'*monos, w22'*monos];
      
  %scaling factor rho is a polynomial of degree 3
 rho_coeffs = sdpvar(length(monos),1); % decision variables
 rho = rho_coeffs'*monos;
 
% %  find dw/dt, assuming partial dw/dt =0
 w11dot = jacobian(W(1,1),[x,y]) * fx;
 w12dot = jacobian(W(1,2),[x,y]) * fx;
 w22dot = jacobian(W(2,2),[x,y]) * fx;
 Wdot = [w11dot, w12dot; 
         w12dot, w22dot];
%  
%  w11 = sdpvar(1,1);
%  w12 = sdpvar(1,1);
%  w22 = sdpvar(1,1);
%  W = [w11,w12;w12,w22];
%  Wdot = zeros(2,2);
 
 alpha1 = 0.1;
 alpha2 = 100;
 lambda = 5;
 
 % there are four sos constraints corresponding to the CCM condition, W >
 % alpha1, W < alpha2, rho > 0
 dummy = sdpvar(length(W),1);
 F = [ sos(-dummy'*(-Wdot+ W*A' + A*W - rho*(B*B') + lambda*W)*dummy);
       sos(dummy' *(W  -  alpha1*eye(size(W)))*dummy);
       sos(dummy' *(alpha2*eye(size(W)) - W)*dummy);
       sos(rho)];
   

 solvesos(F, 0,[],[w11;w12;w22;rho_coeffs]);
   
W2 = replace(W,w11,value(w11));
W2 = replace(W2,w12,value(w12));
W2 = replace(W2,w22,value(w22));
  

%invert
det = (W2(1,1)*W2(2,2) - W2(1,2).^2)
M = [W2(2,2) , -W2(1,2);
      -W2(1,2), W2(1,1)]  ;

% M = eye(size(W2))/W2;

 %%
 % If we assume a flat metric and look at states near the origin 
 % we can calcualte the static feedback gain so that u = Kx
 

x0 = sdpvar(2,1);
sdpvar s;

xstar = [0;0];
ustar = 0;
    
trans = x0*s + (1-s)*xstar;
det = replace(det,x,trans(1));
det = replace(det,y,trans(2));

dgds = x0-xstar;

p = -0.5*rho*B'*M * dgds ; 

%replace symbolic expressions with values in order to speed up integration.
p = replace(p,rho_coeffs,value(rho_coeffs));
p = replace(p,w11,value(w11));
p = replace(p,w12,value(w12));
p = replace(p,w22,value(w22));

q = replace(p, x,trans(1)); %convert from (x,y) to s
q = replace(q, y,trans(2));


dyn = @(xt,ut,t) [-xt(2) - 1.5*xt(1)^2 - 0.5*xt(1)^3; 
            xt(1) + ut];
h = 0.05;
z = [-1; 2];
t = 0;

t= zeros(h*300,1);
for k = 1:300
    
    q2 = replace(q,x0,z(:,k));
    det2 = replace(det,x0,z(:,k));
    ut = 0;
    h = 0.01;
    for j = 1 : 100
        ut = ut + h*replace(q2,s, j*h)/replace(det2,s,j*h);
    end
%     ut = int(q2/det,s,0,1);
    z(:,k+1) = RungeKutta4(dyn, z(:,k), ut, k*h ,h   );
    t(k+1) = t(k)+h;
end
 
    
 plot(t,z(1,:))
 hold on
 plot(t,z(2,:))
 

 
 
 