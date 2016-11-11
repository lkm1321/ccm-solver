function xdot = LVcontrolle_withoutTransvers(~, x)

x1 = x(1);
x2 = x(2); 
u = - (2.2728e+03*x1+1.0708e+03*x2-39.5344*x1*x2-0.1071*x1^2-3.7084e+03*x2^2); 

xdot = [x1-x1*x2 + u;
            -x2+x1*x2]; 

end 