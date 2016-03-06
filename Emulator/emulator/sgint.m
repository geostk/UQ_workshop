function sgvalue = sgint(l,u,w,theta,ctype)
%%%% Compute the integration of correlation function 
%%%% \int_l^u R (x,w)g (x)dx where g (x)~U[l,u]
%%%% Input l         : lower limit 
%%%%       u         : upper limit
%%%%       w         : given constant
%%%%       theta     : estimated scale parameter in the correlation function
%%%%       ctype     : correlation family (0: Gaussian, 2: Cubic)
%%%% Output sgvalue  : answer of the integration


%%% Compute the integration for Gaussian correlation function
if ctype == 0 
    temp = normcdf ((u-w)*sqrt (2*theta)) - normcdf ((l-w)*sqrt (2*theta));
    sgvalue = sqrt (pi/theta) / (u-l) * temp;
    
%%% Compute the integration for Cubic correlation function    
elseif ctype == 2

    lstar = max(w-theta,l);
    ustar = min(w+theta,u);
    
    if l > (w - theta/2) & u <= (w + theta/2)
	sgvalue = (u-l) - 2*((u-w)^3-(l-w)^3)/theta^2 + 3*((u-w)^4 + (l-w)^4)/(2*theta^3);
        
    elseif l > (w - theta/2) & u > (w + theta/2)
	sgvalue = (w-l) + 3*theta/8 + 2*(l-w)^3/theta^2 + (3*(l-w)^4-(theta-ustar+w)^4)/(2*theta^3);
        
    elseif l <= (w - theta/2) & u <= (w + theta/2)
	sgvalue = (u-w) + 3*theta/8 - 2*(u-w)^3/theta^2 + (3*(u-w)^4-(theta+lstar-w)^4)/(2*theta^3);
        
    elseif l <= (w - theta/2) & u > (w + theta/2)
	sgvalue = 3*theta/4 - ((theta-ustar+w)^4 + (theta+lstar-w)^4)/(2*theta^3);
    end
    
    sgvalue = sgvalue / (u-l);
end


% example
% l=0;u=1;
% w=0.6;theta=1.6;
% w=0.2;theta=0.8;
% w=0.9;theta=0.8;
% w=0.7;theta=0.2;

%% numerical integration
% if l > (w - theta/2) & u <= (w + theta/2)
%         inte = quad(@(x) 1-6*((x-w)./theta).^2 + 6*(abs(x-w)./theta).^3, l,u) ;
% 
% elseif l > (w - theta/2) & u > (w + theta/2)
%         inte1 = quad(@(x) 1-6*((x-w)./theta).^2 + 6*(abs(x-w)./theta).^3, l,w+theta/2) ;
%         inte2 = quad(@(x) 2*(1-abs(x-w)./theta).^3, w+theta/2,ustar);
%         inte = inte1+inte2;
%         
% elseif l <= (w - theta/2) & u <= (w + theta/2)
%         inte1 = quad(@(x) 2*(1-abs(x-w)./theta).^3, lstar,w-theta/2);
%         inte2 = quad(@(x) 1-6*((x-w)./theta).^2 + 6*(abs(x-w)./theta).^3, w-theta/2,ustar);
%         inte = inte1+inte2;
%         
% elseif l <= (w - theta/2) & u > (w + theta/2)
%         inte1 = quad(@(x) 2*(1-abs(x-w)./theta).^3, lstar,w-theta/2);
%         inte2 = quad(@(x) 1-6*((x-w)./theta).^2 + 6*(abs(x-w)./theta).^3, w-theta/2,w+theta/2);
%         inte3 = quad(@(x) 2*(1-abs(x-w)./theta).^3, w+theta/2,ustar);
%         inte = inte1+inte2+inte3;
% end
