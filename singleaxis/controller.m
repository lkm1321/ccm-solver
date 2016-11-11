function mu = controller(theta,omega,k)

persistent theta_final
persistent alpha; 
persistent c; 
persistent mu_s; 

if (isempty(theta_final))
    theta_final = -pi/4;     
end

if (isempty(alpha))
   alpha = 0.07;  
end

if (isempty(c))
    c = 0.07; 
end

if (isempty(mu_s))
   mu_s = 0.12; 
end

s = omega + alpha/2*(theta-theta_final-(sin(2*theta)/2-sin(2*theta_final)/2));
mu = -(alpha/k*omega*sin(theta)+c*sign(s)*sign(sin(theta))); 

if (abs(mu)>mu_s)
    mu = sign(mu)*mu_s; 
end

end