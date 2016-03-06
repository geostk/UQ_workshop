function sens = sensitivity(x,y,beta,sigma2,invSigma,theta,ctype);
%%%% Estimate the sensitivity indices
%%%% Input x         : n by p design matrix
%%%%       y         : n by 1 output vector  
%%%%       beta      : estimated mean parameter
%%%%       sigma2    : esitmated process variance
%%%%       invSigma  : n by n inverse covariance function corresponding to 
%%%%                   the training data and the estimated correlation parameters
%%%%       theta     : estimated correlation parameters
%%%%       ctype     : correlation family (0: Gaussian, 2: Cubic)
%%%% Output sens     : a structure of the estimated sensitivity indices
%%%%        sens.sme   : estimated main effect sensitivity indices
%%%%        sens.ste   : estimated total effect sensitivity indices
%%%%        sens.totalvar   : estimated total variance
%%%%        sens.ngrid : number of grids of the inputs
%%%%        sens.maineffect : main effect values evaluated at the grids
%%%% It calls the function
%%%%        sgint.m
%%%%        dbint.m
%%%%        Cfunc.m which calls sgint.m, dbint.m, mxint.m

%%% dimension of design
[n, p] = size(x);

%%% input region
l = min(x);
u = max(x);

%%%%%%%%%%%%%%%%%%% Computation of Sensitivity Indices %%%%%%%%%%%%%%%%%%%%%
%%% n by 1 vector q
q = zeros(n,1);
for i=1:n
    temps = 1;
    for k=1:p
        temps = temps * sgint(l(k),u(k),x(i,k),theta(k),ctype);
    end
    q(i) = sigma2 * temps;
end

%%% The third component of variance is fixed for all nputs, so compute early
tempr = 1 ;
for k=1:p
    tempr = tempr * dbint(l(k),u(k),theta(k),ctype);
end
Vs3 = sigma2 * tempr - trace(invSigma*q*q');

%%%%%%%%%%% Compute total variance %%%%%%%%%%%%%
inc = 1:p;
exc = [];

%%% matrix C 
C = Cfunc(n,l,u,x,inc,exc,theta,sigma2,ctype);

%%% Vs1
Vs1 = sigma2 - trace(invSigma * C);

%%% Vs2
Vs2 = (y-beta)' * invSigma * (C - q*q') * invSigma * (y-beta);

%%% total variance
vt = Vs1 + Vs2 - Vs3;

%%%%%% Compute the main effect variance and total effect variance %%%%%
all = 1:p;
Vs = zeros(2,p);
for z = 1:p
    for index = 1:2
        if index == 1                % for sme (main effect)
             inc = z;                  
             exc = all (all ~= inc);   
         else                        % for ste (total effect)
             exc = z;                  
             inc = all (all ~= exc);   
         end        
    
         %%% matrix C
         C = Cfunc(n,l,u,x,inc,exc,theta,sigma2,ctype);
    
         %%% Vs1
         temp1 = 1;
         for kk=1:length(exc)
             k = exc(kk);
             temp1 = temp1 * dbint(l(k),u(k),theta(k),ctype);
         end
         Vs1 = sigma2 * temp1 - trace(invSigma * C);

         %%% Vs2
         Vs2 = (y-beta)' * invSigma * (C - q*q') * invSigma * (y-beta);

         %%% Vs
         Vs(index,z) = Vs1 + Vs2 - Vs3;
     end    
end

%%% main effect variance
vsme = Vs(1,:);
%%% main effect sensitivity index
sme = vsme / vt;

%%% total effect variance
vste = vt - Vs(2,:);
%%% total effect sensitivity index
ste = vste / vt;


%%%%%%%%%%%%%%%%%%% Computation of Main Effect Function %%%%%%%%%%%%%%%%%%%%%
%%% make grids for each input
ngrid = 21;
xe = zeros(ngrid,p);
for i = 1:p
    xe(:,i) = [l(i) : (u(i)-l(i))/(ngrid-1) : u(i)]';
end

%%% evaluate the main effect function based on the grid values
all = 1:p;
covgrid = zeros(ngrid,n);
effect = zeros(ngrid,p);

for z = 1:p
    
    inc = z;
    exc = all (all ~= inc);   
    
    for w = 1:ngrid
        for i = 1:n
            
            tempca = 1;
            for kk = 1:length(exc)
                k = exc(kk);
                tempca = tempca * sgint(l(k),u(k),x(i,k),theta(k),ctype);
            end
            
            if ctype == 0
                tempcb = exp(-theta(z)*(xe(w,z)-x(i,z)).^2);
            elseif ctype == 2
                dist = abs(xe(w,z)-x(i,z));
                if dist <= (theta(z)/2)
                    tempcb = 1-6*(dist/theta(z))^2 + 6*(dist/theta(z))^3;
                elseif dist > (theta(z)/2) & dist <= theta(z)
                    tempcb = 2*(1-dist/theta(z))^3;
                elseif dist > theta(z)
                    tempcb = 0;
                end
            end
            
            covgrid(w,i) = sigma2 * tempca * tempcb;
        end
        
        effect(w,z) = beta + covgrid(w,:)*invSigma*(y-beta);    
    end    
end

%%% normalized grids of the inputs
xgrid = [0 : 1/(ngrid-1) : 1]';

labels = strcat('x_', num2str([1:p]'));
plot(xgrid, effect); 
xlabel('Scaled Input'); ylabel('Main Effect');
title ('Main Effect Plots');
legend(labels, 'Location', 'Best');


print mainplot '-djpeg90';close;

%%% save the results
sens.sme = sme;
sens.ste = ste;
sens.totalvar = vt;
sens.ngrid = ngrid;
sens.maineffect = effect;


  