function [loglik] = llik1d(corparms, Y, X, ctype, fittype, ...
                       F, dimcp);
%%%% Compute the likelihood function
%%%% Input:     'corparms' == the dimcp vector of correlation parameters
%%%%            'Y'        == n by 1 vector of responses 
%%%%            'X'        == n by p matrix of inputs for the Y vector
%%%%            'ctype'    == the correlation type: 0==Gaussian Correlation, 1==Prod.Exp., 
%%%%            'fittype'  == 0 if ML and 1 if REML
%%%%            'F'        == the F matrix for the n observations
%%%%            'dimcp'    == correlation parameter vector
%%%% Output:    negative likelihood corresponding to the inputs. 

n = length(Y);
nvar = length(X(1,:));
flag = 0;
npar = length(F(1,:));
p=dimcp/2;

% Calculate the correlation matrix of Y  
if ctype   ==   2,
    %% cubic correlation matrix: anisotropic correlation function
    [n,p]   =   size(X);
    R       =   zeros(n,n);
    for i   =   1:n
        for j        =   1:n
            R(i,j)   =   cubiccorrfn(X(i,:)-X(j,:),corparms(1:p));
        end
    end
else 
    %% Gaussian or PowerExponential function
    R = cormatexp(X,corparms(1:p),corparms(p+1:2*p));
end
[cholR, p] = chol(R);
R = R;% * 0.99 + 0.01 * eye(n,n); <-- Use the commented code to fit a model
      %                               with random error.

% Compute the likelihood or restricted likelihood value    
if (rcond(R) > 1e-12 & p==0)  %if not this then bad corparms
    RinvF = R\F;
    tFRinvF = F'*RinvF;
    RinvY = R\Y;
    tFRinvFinvFt = tFRinvF\F';
    Betahat = tFRinvFinvFt*RinvY;
    resid = Y - F*Betahat;
    Rinvresid = R\resid;
    detR = det(R);
    if (detR <= 0.0)
        loglik = 1.0e10;
    else
        if (fittype == 0)  
            T1hat2 = (resid'*Rinvresid)/n;
            loglik = n*log(T1hat2) + log(detR) + n + n*log(2*pi);
            %C-PERK used n*log(T1hat2) + log(detR). 
        elseif (fittype == 1)
            T1hat2 = (resid'*Rinvresid)/(n - npar); 
            %C-PERK used (n-npar-2).
            detNmat = det(tFRinvF);
            loglik = (n-npar)*log(T1hat2) + log(detNmat) + log(detR);
            %C-PERK used (n-npar)*log(T1hat2) + log(detNmat) + log(detR) 
            npar;
            a1=(n-npar)*log(T1hat2); 
            b1=log(detNmat);
            c1=log(detR);
        end
        if (T1hat2 < 0.0)
            loglik = 1.0e10;
        end
    end
else
    loglik = 1.0e10;
end    
