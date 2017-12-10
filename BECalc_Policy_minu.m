% Function used to find the Bellman error, given the value function estimate, the reward, discount factor, number of state action pairs, and the noise information (this essentially gives us the information on the transition probabilty matrix)

function BE = BECalc_Policy_minu(Valfun,R,beta,StaLen,noise)

BE = zeros(1,StaLen);


for x = 1:StaLen

    
    RplusValfun_x = zeros(1,StaLen);
    
    Valfun_xtp1 = 0;
    NumOfXtp1 = 0;
    
    % Look at all possible xtp1, give xt = x, to calculate the h(xtp1) in the case that noise "happened"
    
    for xtp1 = 1:StaLen
        if(R(x,xtp1) ~= -100)

        Valfun_xtp1 = Valfun_xtp1 + Valfun(xtp1);
        NumOfXtp1 = NumOfXtp1 + 1;
        
        end
    end
    
    
    E_Valfun_xtp1 = Valfun_xtp1./NumOfXtp1;
    
    
    for a = 1:StaLen
        if(R(x,a) ~= -100)
            
            xdash = a;
            RplusValfun_x(a) = R(x,a) + beta*((1-noise)*Valfun(xdash)  +  noise*E_Valfun_xtp1);       
            
        end
    end
    
    
    BE(x) = max(RplusValfun_x) - Valfun(x);

end    
