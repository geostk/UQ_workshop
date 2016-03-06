function [corparms, loglik] = estgaussiancor(x,y,F,corparms,initnum);
%%%% Estimate the correlation parameters in the Gaussian Correlation Function 
%%%% Input x         : n by p matrix x, 
%%%%       y         : n by 1 matrix y
%%%%       corparms  : a structure of the correlation parameters 
%%%%       corparms.scale: p by 1 parameters in the cubic correlation
%%%%       function
%%%% Output corparms : estimated correlation parameters
%%%%        loglik   : minimized negative log likelihood.

%%% To caluculate the mindist and maxdist
%%% given x
[n p] = size(x);
selectx = zeros(n*(n-1)/2,p);
%%% save the distances in the selectx matrix
count = 1;
for i=1:n-1,
   for j=(i+1):n,
      selectx(count,:) = (x(i,:)-x(j,:)).^2;
      count = count + 1;
   end
end
% set a bound for the minimum distance
selectx(selectx < 0.33/n) = 0.33/n;
        
minx = min(selectx);
maxx  = max(selectx);
lbd  = zeros(p,1);
ubd  = zeros(p,1);

for i=1:p,
    ubd(i) = -log(0.01^(1/p))/minx(i);
    lbd(i) = -log(0.99^(1/p))/maxx(i);
end

%%% Get lower bound and upper bound
lowbd   =   lbd;
highbd  =   ubd;

%%% Make the LHS design for the scale parameters %%%
scalepar=lhsdesign(200*p,p);
for i = 1:p,
    scalepar(:,i) = (highbd(i)-lowbd(i)) * scalepar(:,i) ...
        + lowbd(i);
end  
corend=zeros(200*p,p);
loglikarray=zeros(200*p,1);


%%% Calculate likelihood/restricted likelihood for all the initial designs
if corparms.fittype==0, % MLE
    for i=1:200*p,
      loglikarray(i) = calicatMLE(scalepar(i,1:p),x,y,F,0);
    end
elseif corparms.fittype==1, % REML
    for i=1:200*p,
      loglikarray(i) = calicatREML(scalepar(i,1:p),x,y,F,0);
    end
end


%%% Determine the best initnum*p starting values
% Sort loglikarray in ascending order
[sortedLogLik index] = sort(loglikarray);
% Choose the smallest values 
loglikarray = sortedLogLik(1:initnum*p, :);
% Choose the correlation parameter values which produced
% the samllest values of loglik function. These values will be 
% used as starting values for fmincon optimization.
scalepar = scalepar(index(1:initnum*p, :), :);

%%% Let the search be in a wider range.
lowbd(1:p)=lbd/50;
highbd(1:p)=ubd*50;

% Set the options for the fmincon function.
options = optimset('Display','off','LargeScale','on',...
                   'TolFun',.0001, 'Algorithm','active-set');

%%% Check the estimation: MLE or REML %%% 
if corparms.fittype==0, %%% MLE,
    for i=1:initnum*p,
      corparinit=scalepar(i,:);
      [cort,fvals] = fmincon(@(cort)calicatMLE(cort,x,y,F,0),corparinit,...
                      [],[],[],[],lowbd,highbd,[],options);
      corend(i,:) = cort;
      loglikarray(i) = fvals;
    end 
elseif corparms.fittype==1,
    for i=1:initnum*p,
        corparinit=scalepar(i,:);
        [cort,fvals] = fmincon(@(cort)calicatREML(cort,x,y,F,0),...
                        corparinit,...
                        [],[],[],[],lowbd,highbd,[],options);
        corend(i,:) = cort;
        loglikarray(i) = fvals;
    end
end

%%% Find the best answer and save output
%%% As we minimize the negative likelihood, we should get the minimum.
[loglik index]= min(loglikarray);
corparms.scale=corend(index,:);
corparms.smoothness=2*ones(p,1);

