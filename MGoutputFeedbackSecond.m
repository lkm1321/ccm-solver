function xdot = MGoutputFeedbackSecond(~, x, xstar1)

    xhat1 = x(1); 
    xhat2 = x(2); 
    
    y = x(4) + 0.2*randn(1); 
    
    feedBackTop = (364.3177*xhat1+35.4547*y+261.2955*xhat1^3+26.3638*xhat1^2+2.9268*xhat1*y-0.0919*y^2+59.3999*xhat1^2*y+4.5435*xhat1*y^2+0.1212*y^3); 
    feedBackBottom = -(364.3177*xhat1+35.4547*xhat2+261.2955*xhat1^3+26.3638*xhat1^2+2.9268*xhat1*xhat2-0.0919*xhat2^2+59.3999*xhat1^2*xhat2+4.5435*xhat1*xhat2^2+0.1212*xhat2^3); 

    ctrl = MGgetControl(x(1:2), xstar1); 
    xhatdot = MGopenloop(0, x(1:2), ctrl) + feedBackTop + feedBackBottom; 
    xtruedot = MGopenloop(0, x(3:4), ctrl); 
    xdot = [xhatdot; xtruedot]; 
    
end