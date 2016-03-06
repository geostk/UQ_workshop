function Cmatrix = Cfunc(n,l,u,x,inc,exc,theta,sigma2,ctype)
%%%% Compute \int_l^u Cov_p[\eta_s(x_s),y_sim]' Cov_p[\eta_s(x_s),y_sim] g(x_s)dx_s
%%%% Input l         : lower limit 
%%%%       u         : upper limit
%%%%       x         : n by p design matrix
%%%%       inc       : set of inputs for which compute sensitivity indices
%%%%       exc       : complement of inc
%%%%       theta     : estimate scale parameter in the correlation function
%%%%       sigma2    : estimated process variance
%%%%       ctype     : correlation family (0: Gaussian, 2: Cubic)
%%%% Output Cmatrix  : n by n matrix

%%% Compute the elements of the matrix C by calling the function sgint and mxint 
Cmatrix = zeros(n,n);
for i=1:n
    for j=1:n            
        tempa = 1;
        if isempty(exc)==0
            for kk=1:length(exc)
                k = exc(kk);
                tempa = tempa * sgint(l(k),u(k),x(i,k),theta(k),ctype) * sgint(l(k),u(k),x(j,k),theta(k),ctype);
            end
        end    

        tempb = 1;
        if isempty(inc)==0
            for kk=1:length(inc)
                k = inc(kk);
                tempb = tempb * mxint(l(k),u(k),x(i,k),x(j,k),theta(k),ctype);
             end    
        end
         
        Cmatrix(i,j) = sigma2^2 * tempa * tempb;
    end    
end