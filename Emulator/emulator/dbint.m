function dbvalue = dbint_mod(l,u,theta,ctype)
%%%% Compute the double integration of the correlation function
%%%% \int_l^u \int_l^u R (x1,x2)g (x1)g (x2)dx1dx2 where g (x1) and g (x2) ~ U[l,u]
%%%% Input l         : lower limit
%%%%       u         : upper limit
%%%%       theta     : estimated scale parameter in the correlation function
%%%%       ctype     : correlation family (0: Gaussian, 2: Cubic)
%%%% Output dbvalue  : answer of the double integration

%%% Compute the double integration for Gaussian correlation function
if ctype == 0
    
    diff = u-l;
    temp1 = (1/theta) * ( sqrt (2*pi) * normpdf (diff*sqrt (2*theta)) -1 );
    temp2 = diff * sqrt (pi/theta) * (2*normcdf (diff*sqrt (2*theta)) -1 );
    dbvalue = (temp1 + temp2)/diff^2;

%%% Compute the double integration for Cubic correlation function
elseif ctype == 2 
    
    if (u-l) > theta
        dbvalue = (23*theta^2)/40 + 3/4*theta*(-l - theta + u);
        
    elseif (u-l) <= theta & (u-l) > (theta/2)
        dbvalue = (l+theta-u)^5/(5*theta^3) +(3/4)*theta*(u-l) -(7/40)*theta^2;
        
    elseif (u-l) <= (theta/2)
        dbvalue = (l-u)^2 - (l - u)^4*(3*l - 3*u + 5*theta)/(5*theta^3);
    end
    
    dbvalue = dbvalue / (u-l)^2;
end    



% scaled input
% if (u-l) > theta
%     dbvalue = theta^2*23/40;
% elseif (u-l) <= theta & (u-l) > (theta/2)
%     dbvalue = 2 - 1/(5*theta^3) + 1/theta^2 - 2/theta - theta/4 + theta^2/40;
% elseif (u-l) <= (theta/2)
%     dbvalue = 1 - 1/theta^2 + 3/(5*theta^3);
% end

