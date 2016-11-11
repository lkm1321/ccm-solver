function xdot = MGopenloop(~, xcurr, ustar)

    x1 = xcurr(1);
    x2 = xcurr(2); 
    xdot =[-0.5*x1^3-1.5*x1^2-x2;... 
           x1+ustar]; 

end