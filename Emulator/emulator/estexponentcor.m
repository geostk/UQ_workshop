function [corparms,loglik] = estexponentcor(x,y,F,fittype,initnum);
%%%% Estimate the correlation parameters in the Power Exponential 
%%%%       Correlation Function 
%%%% Inputs 
%%%%       x         : n by p matrix of inputs 
%%%%       y         : n by 1 column vector of outputs
%%%%       fittype   : (=0 for MLE and =1 for RML),  
%%%%       initnum:  : numbet \times p of starting points used in
%%%%                        maximization of the loglikelihood
%%%% Outputs 
%%%%        corparms: a structure of the correlation parameters with 
%%%%                    structure elements
%%%%           corparms.scale = estimated thetas
%%%%           corparms.smoothness = estimated powers
%%%%       loglik   : minimized negative log likelihood.

%%% Caluculate the minimum and maximum Euclidean 
%%%  interpoint distance among rows inputs of  x
[n p] = size(x);
selectx = zeros(nchoosek(n,2),p);
%%% save the distances in the selectx matrix
count = 1;
for i=1:n-1,
   for j=(i+1):n,
      selectx(count,:) = (x(i,:)-x(j,:)).^2;
      count = count + 1;
   end
end
% Set a lower bound for the minimum distance
selectx(selectx < 0.33/n) = 0.33/n;

minx = min(selectx);
maxx  = max(selectx);
lbd  = zeros(p,1);
ubd  = zeros(p,1);

% Construct lower and upper bounds for each \theta
for i=1:p,
    ubd(i) = -log(0.01^(1/p))/minx(i);
    lbd(i) = -log(0.99^(1/p))/maxx(i);
end

% Make the LHS design for the scale parameters %%%
lowbd=zeros(2*p,1);
highbd=zeros(2*p,1);
lowbd(1:p)=lbd;
highbd(1:p)=ubd;
highbd(p+1:2*p)=2;

max_p = 20*p;
scalepar=lhsdesign(max_p,2*p,'iterations',50);
for i = 1:2*p,
  scalepar(:,i) = (highbd(i)-lowbd(i)) * scalepar(:,i) ...
                   + lowbd(i);
end  
corend=zeros(max_p,2*p);
loglikarray=zeros(max_p,1);

% Calculate likelihood/restricted likelihood 
%   for all max_p candidate designs %%% 
if fittype==0, %%% MLE,
    for i=1:max_p,
        loglikarray(i) = calicatMLE(scalepar(i,:),x,y,F,1);
    end
elseif fittype==1,
    for i=1:max_p,
        loglikarray(i) = calicatREML(scalepar(i,:),x,y,F,1);
    end
end  % end if corparms.fittype==0

% Determine the best initnum*p starting values
%   and sort loglikarray in ascending order
[sortedLogLik index] = sort(loglikarray);

% Choose the smallest loglik values 
loglikarray = sortedLogLik(1:initnum*p, :);

% Select the correlation parameter values that produced
%   the samllest values of loglik function. These values will be 
%   used as starting values for fmincon optimization.
scalepar = scalepar(index(1:initnum*p, :), :);

% Widen the \theta search 
lowbd(1:p)=lbd/50;
highbd(1:p)=ubd*50;

% Set the options for the fmincon function.
options = optimset('Display','off','LargeScale','on',...
                   'TolFun',.0001, 'Algorithm','active-set');
         
% Check the estimation: MLE or REML %%% 
if fittype==0, %%% MLE,
    for i=1:initnum*p,
      corparinit=scalepar(i,:);
      [cort,fvals] = fmincon(@(cort)calicatMLE(cort,x,y,F,1),...
            corparinit,[],[],[],[],lowbd,highbd,[],options);
      corend(i,:) = cort;
      loglikarray(i) = fvals; 
    end    
elseif fittype==1,
    for i=1:initnum*p,
        corparinit=scalepar(i,:);
        [cort,fvals] = fmincon(@(cort)calicatREML(cort,x,y,F,1),...
            corparinit,[],[],[],[],lowbd,highbd,[],options);
        corend(i,:)     = cort;
        loglikarray(i)  = fvals; 
    end
end

% Find the set of correlation and power parameters that
%     minimize the negative likelihood answer and save output
[loglik index]= min(loglikarray);
corparms.scale = corend(index,1:p);
corparms.smoothness = corend(index,p+1:2*p);
end % function
