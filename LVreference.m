function xdot = LVreference(~, xcurr, xstar2)

x1 = xcurr(1);
x2 = xcurr(2); 
xstar1 = 1; 

utop = - (740.0241*x1-212.5183*x2+22.1967*x1*x2-5.4155*x1^2+894.4307*x2^2-0.1447*x1^2*x2+112.1987*x1*x2^2-1.5439*x2^3); 
ubottom = (740.0241*xstar1-212.5183*xstar2+22.1967*xstar1*xstar2-5.4155*xstar1^2+894.4307*xstar2^2-0.1447*xstar1^2*xstar2+112.1987*xstar1*xstar2^2-1.5439*xstar2^3);
ustar = -xstar1 +xstar1*xstar2; 

% u = 1; 
xdot = [x1-x1*x2 + utop + ubottom + ustar;
            -x2+x1*x2]; 



end 