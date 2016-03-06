function [corparms, loglik] = estcubiccor(x,y,F,corparms,initnum);
%%%% Estimate the correlation parameters for the cubic correlation function
%%%% using either MLE or REML
%%%% Input x         : n by p matrix x, 
%%%%       y         : n by 1 matrix y
%%%%       corparms  : a structure of the correlation parameters 
%%%%       corparms.scale: p by 1 parameters in the cubic correlation function
%%%%       
%%%% Output corparms : estimated correlation parameters
%%%%        loglik   : minimized negative log likelihood.


%%% Make the LHS design for the scale parameters %%%
p=size(x,2);
scalepar=lhsdesign(200*p,p);
scalepar=(100/p)*scalepar;
corend=zeros(initnum*p,p);
loglikarray=zeros(200*p,1);
lowbd=zeros(p,1);
highbd=lowbd;
lowbd(:)=0.00001;
highbd(:)=100/p;


%%% Calculate likelihood/restricted likelihood for all the initial designs %%% 
if corparms.fittype==0, %%% MLE,
    for i=1:200*p,
      loglikarray(i) = calicatMLE(scalepar(i,1:p),x,y,F,2);
    end
elseif corparms.fittype==1,
    for i=1:200*p,
      loglikarray(i) = calicatREML(scalepar(i,1:p),x,y,F,2);
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

% Set the options for the fmincon function.
options = optimset('Display','off','LargeScale','on',...
                   'TolFun',.0001, 'Algorithm','active-set');
               
%%% Check the estimation: MLE or REML %%% 
if corparms.fittype==0, %%% MLE,
    for i=1:initnum*p,
      corparinit=scalepar(i,1:p);
      [cort,fvals] = fmincon(@(cort)calicatMLE(cort,x,y,F,2),corparinit,...
                      [],[],[],[],lowbd/50,highbd*500,[],options);
      corend(i,:) = cort;
      loglikarray(i) = fvals;
    end
elseif corparms.fittype==1,
    for i=1:initnum*p,
        corparinit=scalepar(i,1:p);
        [cort,fvals] = fmincon...
            (@(cort)calicatREML(cort,x,y,F,2),corparinit,...
                        [],[],[],[],lowbd/50,highbd*500,[],options);
        corend(i,:) = cort;
        loglikarray(i) = fvals;
    end
end

%%% find the best answer and save output
[loglik index]= min(loglikarray);
corparms.scale=corend(index,:);



