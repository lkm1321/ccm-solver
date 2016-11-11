function xdot = LVModel(~, xcurr)

x1 = xcurr(1);
x2 = xcurr(2); 
u = - (740.0241*x1-212.5183*x2+22.1967*x1*x2-5.4155*x1^2+894.4307*x2^2-0.1447*x1^2*x2+112.1987*x1*x2^2-1.5439*x2^3); 

xdot = [x1-x1*x2 + u;
            -x2+x1*x2]; 
end 