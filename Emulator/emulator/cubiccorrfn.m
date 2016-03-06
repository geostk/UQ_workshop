function cubiccorr = cubiccorrfn(x,theta);
%%%% Computes the value of the cubic correlation function
%%%% given parameter theta (a px1 vector) at x (a 1xp point)

p=length(theta);
r=1;
for i=1:p
    if abs(x(i))>abs(theta(i)) 
        c=0; 
    elseif abs(x(i))>(abs(theta(i))/2)
        c=2*(1-(abs(x(i))/abs(theta(i))))^3; 
    else
        c=1-6*(x(i)/abs(theta(i)))^2 + 6*((abs(x(i))/abs(theta(i)))^3); 
    end;
    r=r*c;
end;
cubiccorr=r;