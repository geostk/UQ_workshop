function mxvalue = mxint_mod(l,u,w1,w2,theta,ctype)
%%%% Compute the integration of the product of the correlation functions
%%%% \int_l^u R(x,w1)R(x,w2)g(x)dx where g(x)~U[l,u]
%%%% Input l         : lower limit 
%%%%       u         : upper limit
%%%%       w1        : first given constant
%%%%       w2        : second given constant
%%%%       theta     : estimated scale parameter in the correlation function
%%%%       ctype     : correlation family (0: Gaussian, 2: Cubic)
%%%% Output mxvalue  : answer of the integration

%%% Compute the integration for Gaussian correlation function
if ctype == 0
    tempa1 = exp(-(1/2)*theta*(w1-w2)^2);
    tempa2 = sgint(l,u,(w1+w2)/2,2*theta,ctype);
    mxvalue = tempa1 * tempa2;
    
%%% Compute the integration for Cubic correlation function
elseif ctype == 2

    delta = abs(w1-w2);
    w1=min(w1,w2);
    w2=w1+delta;
    lstar = max(w2-theta,l);
    ustar = min(w1+theta,u);
    
    %numberical integration
    if delta == 0 
        if l > (w1-theta/2)
            if u <= (w1+theta/2)
                inte = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                 (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), l,u);
            elseif u > (w1+theta/2) 
                inte1 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), l,w1+theta/2);
                inte2 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                                  (2*(1-abs(x-w2)./theta).^3), w1+theta/2,ustar);
                inte = inte1+inte2;              
            end
        elseif l <= (w1-theta/2)    
             if u <= (w1+theta/2)
                inte1 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                                  (2*(1-abs(x-w2)./theta).^3), lstar,w1-theta/2);
                inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w1-theta/2,u);
                inte = inte1+inte2;              
            elseif u > (w1+theta/2)
                inte1 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                                  (2*(1-abs(x-w2)./theta).^3), lstar,w1-theta/2);
                inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w1-theta/2,w1+theta/2);
                inte3 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                                  (2*(1-abs(x-w2)./theta).^3), w1+theta/2,ustar);
                inte = inte1+inte2+inte3;              
            end
        end   
    elseif (delta > 0) & (delta <= (theta/2)) 
        if l > (w2-theta/2) 
            if u <= (w1+theta/2)
                inte = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                 (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), l,u);
            elseif (u > (w1+theta/2)) & (u <= (w2+theta/2))
                inte1 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), l,w1+theta/2);
                inte2 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w1+theta/2,u);
                inte = inte1+inte2;              
            elseif u > (w2+theta/2) 
                inte1 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), l,w1+theta/2);
                inte2 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w1+theta/2,w2+theta/2) ;
                inte3 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                                  (2*(1-abs(x-w2)./theta).^3), w2+theta/2,ustar);
                inte = inte1+inte2+inte3;
            end
        elseif l <= (w2-theta/2) & l > (w1-theta/2)
            if u <= (w1+theta/2)
                inte1 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (2*(1-abs(x-w2)./theta).^3), l,w2-theta/2);
                inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,u);
                inte = inte1+inte2;              
            elseif (u > (w1+theta/2)) & (u <= (w2+theta/2))
                inte1 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (2*(1-abs(x-w2)./theta).^3), l,w2-theta/2);
                inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,w1+theta/2);
                inte3 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w1+theta/2,u);
                inte = inte1+inte2+inte3;              
            elseif u > (w2+theta/2) 
                inte1 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (2*(1-abs(x-w2)./theta).^3), l,w2-theta/2);
                inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,w1+theta/2);
                inte3 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w1+theta/2,w2+theta/2);
                inte4 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                                  (2*(1-abs(x-w2)./theta).^3), w2+theta/2,ustar) ;
                inte = inte1+inte2+inte3+inte4;
            end
        elseif l < (w1-theta/2)
            if u <= (w1+theta/2)
                inte1 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                                  (2*(1-abs(x-w2)./theta).^3), lstar,w1-theta/2);          
                inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (2*(1-abs(x-w2)./theta).^3), w1-theta/2,w2-theta/2);
                inte3 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,u);
                inte = inte1+inte2+inte3;              
            elseif (u > (w1+theta/2)) & (u <= (w2+theta/2))
                inte1 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                                  (2*(1-abs(x-w2)./theta).^3), lstar,w1-theta/2);          
                inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (2*(1-abs(x-w2)./theta).^3), w1-theta/2,w2-theta/2);
                inte3 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,w1+theta/2);
                inte4 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w1+theta/2,u);
                inte = inte1+inte2+inte3+inte4;              
            elseif u > (w2+theta/2) 
                inte1 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                                  (2*(1-abs(x-w2)./theta).^3), lstar,w1-theta/2) ;          
                inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (2*(1-abs(x-w2)./theta).^3), w1-theta/2,w2-theta/2);
                inte3 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,w1+theta/2);
                inte4 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                                  (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w1+theta/2,w2+theta/2);
                inte5 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                                  (2*(1-abs(x-w2)./theta).^3), w2+theta/2,ustar);
                inte = inte1+inte2+inte3+inte4+inte5;
            end
        end    
    elseif (delta > (theta/2)) & (delta <= theta)
        inte1 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                          (2*(1-abs(x-w2)./theta).^3), lstar,w2-theta/2) ;
        inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                          (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,w1+theta/2);
        inte3 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                          (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w1+theta/2,ustar);
        inte = inte1+inte2+inte3;
    elseif (delta > theta) & (delta <= (3/2*theta))
        inte1 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
                          (2*(1-abs(x-w2)./theta).^3), w2-theta,w1+theta/2);
        inte2 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                          (2*(1-abs(x-w2)./theta).^3), w1+theta/2,w2-theta/2);
        inte3 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
                          (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,w1+theta);
        inte = inte1+inte2+inte3;    
    elseif (delta > (3/2*theta)) & (delta <= (2*theta))
        inte = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*(2*(1-abs(x-w2)./theta).^3), w2-theta,w1+theta);
    elseif delta > (2*theta)
        inte = 0;
    end
    mxvalue = inte / (u-l);
end

        
% example
% l=0;u=3;
% w1=1.1; w2=1.2; theta= 4;
% w1=2.1; w2=1.2; theta= 4;
% w1=2.1; w2=2.2; theta= 4;
% w1=1.1; w2=2.2; theta= 1;
% w1=1.1; w2=3.5; theta= 1.5;


% if delta <= (theta/2) 
%     if l > (w2-theta/2) 
%         if u <= (w1+theta/2)
%             inte = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                              (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), l,u);
%         elseif u > (w1+theta/2) & u <= (w2+theta/2)
%             inte1 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), l,w1+theta/2);
%             inte2 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w1+theta/2,u);
%             inte = inte1+inte2;              
%         elseif u > (w2+theta/2)
%             inte1 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), l,w1+theta/2);
%             inte2 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w1+theta/2,w2+theta/2) ;
%             inte3 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
%                               (2*(1-abs(x-w2)./theta).^3), w2+theta/2,ustar);
%             inte = inte1+inte2+inte3;
%         end
%     elseif l <= (w2-theta/2) & l > (w1-theta/2)
%         if u <= (w1+theta/2)
%             inte1 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (2*(1-abs(x-w2)./theta).^3), l,w2-theta/2);
%             inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,u);
%             inte = inte1+inte2;              
%         elseif u > (w1+theta/2) & u <= (w2+theta/2)
%             inte1 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (2*(1-abs(x-w2)./theta).^3), l,w2-theta/2);
%             inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,w1+theta/2);
%             inte3 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w1+theta/2,u);
%             inte = inte1+inte2+inte3;              
%         elseif u > (w2+theta/2) 
%             inte1 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (2*(1-abs(x-w2)./theta).^3), l,w2-theta/2);
%             inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,w1+theta/2);
%             inte3 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w1+theta/2,w2+theta/2);
%             inte4 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
%                               (2*(1-abs(x-w2)./theta).^3), w2+theta/2,ustar) ;
%             inte = inte1+inte2+inte3+inte4;
%         end
%     elseif l < (w1-theta/2)
%         if u <= (w1+theta/2)
%             inte1 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
%                               (2*(1-abs(x-w2)./theta).^3), lstar,w1-theta/2);          
%             inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (2*(1-abs(x-w2)./theta).^3), w1-theta/2,w2-theta/2);
%             inte3 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,u);
%             inte = inte1+inte2+inte3;              
%         elseif u > (w1+theta/2) & u <= (w2+theta/2)
%             inte1 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
%                               (2*(1-abs(x-w2)./theta).^3), lstar,w1-theta/2);          
%             inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (2*(1-abs(x-w2)./theta).^3), w1-theta/2,w2-theta/2);
%             inte3 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,w1+theta/2);
%             inte4 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w1+theta/2,u);
%             inte = inte1+inte2+inte3+inte4;              
%         elseif u > (w2+theta/2) 
%             inte1 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
%                               (2*(1-abs(x-w2)./theta).^3), lstar,w1-theta/2) ;          
%             inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (2*(1-abs(x-w2)./theta).^3), w1-theta/2,w2-theta/2);
%             inte3 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,w1+theta/2);
%             inte4 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w1+theta/2,w2+theta/2);
%             inte5 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
%                               (2*(1-abs(x-w2)./theta).^3), w2+theta/2,ustar);
%             inte = inte1+inte2+inte3+inte4+inte5;
%         end
%     end    
% elseif delta > (theta/2) & delta <= theta
%         if u <= (w2-theta/2)
%             inte = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                              (2*(1-abs(x-w2)./theta).^3), lstar,u);
%         elseif u > (w2-theta/2) & u <= (w1+theta/2)
%             inte1 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (2*(1-abs(x-w2)./theta).^3), lstar,w2-theta/2);
%             inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,u);
%             inte = inte1+inte2;              
%         elseif u > (w1+theta/2)
%             inte1 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (2*(1-abs(x-w2)./theta).^3), lstar,w2-theta/2) ;
%             inte2 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,w1+theta/2);
%             inte3 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w1+theta/2,ustar);
%             inte = inte1+inte2+inte3;
%         end
% elseif delta > theta & delta <= (3/2*theta)
%             inte1 = quad(@(x) (1-6*((x-w1)./theta).^2 + 6*(abs(x-w1)./theta).^3).*...
%                               (2*(1-abs(x-w2)./theta).^3), w2-theta,w1+theta/2);
%             inte2 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
%                               (2*(1-abs(x-w2)./theta).^3), w1+theta/2,w2-theta/2);
%             inte3 = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*...
%                               (1-6*((x-w2)./theta).^2 + 6*(abs(x-w2)./theta).^3), w2-theta/2,ustar);
%             inte = inte1+inte2+inte3;    
% %           mxvalue = -((1759*theta^7 + 6629*theta^6*(w1 - w2) + 10416*theta^5*(w1 - w2)^2 + 8890*theta^4*(w1 - w2)^3 +...
% %                        4480*theta^3*(w1 - w2)^4 + 1344*theta^2*(w1 - w2)^5 + 224*theta*(w1 - w2)^6 + 16*(w1 - w2)^7)/(560*theta^6)) +...
% %                      ((20*theta^3 - 9*theta^2*(w1 - w2) - 20*theta*(w1 - w2)^2 - 6*(w1 - w2)^3)*(3*theta + 2*w1 - 2*w2)^4)/(560*theta^6);
% elseif delta > (3/2*theta) & delta <= (2*theta)
%              inte = quad(@(x) (2*(1-abs(x-w1)./theta).^3).*(2*(1-abs(x-w2)./theta).^3), w2-theta,ustar);
% %            mxvalue = (2*theta + w1 - w2)^7/(35*theta^6)             
% elseif delta > (2*theta)
%             inte = 0;
% end    