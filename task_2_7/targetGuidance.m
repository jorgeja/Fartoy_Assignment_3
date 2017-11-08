function nu_d  = targetGuidance(nu_t,p_t,tailDistance, U_aMax,delta_pTilde,p)

p_t_tilde = p - p_t; 
pTilde = p_t_tilde - tailDistance*nu_t; %Setpoint is set at a distance tailDistance behind the target. 

 
if sqrt(pTilde'*pTilde) <= 0
    
    nu_d = sqrt(nu_t'*nu_t)*(p_t_tilde/sqrt(p_t_tilde'*p_t_tilde)); % pure pursuit at target speed
    
    if sqrt(p_t_tilde'*p_t_tilde) < tailDistance/2
        nu_d = nu_d/2;
    end
    
else
    kappa = U_aMax*(sqrt(pTilde'*pTilde))/(sqrt(pTilde'*pTilde + delta_pTilde^2));
    nu_a = -kappa*pTilde/(sqrt(pTilde'*pTilde));
    nu_d = nu_a+nu_t;

end

end