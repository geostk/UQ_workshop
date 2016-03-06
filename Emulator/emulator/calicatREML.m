function loglik = calicatREML(corr,X,Y,F,type);
%%%% Calculate the best correlation parameter and the corresponding
%%%% loglikelihood using the logliklihood function and fmincon.

%%%% Set parameters
ctype=type;  % type=0, we use Gaussian Correlation;
             % type=1, we use Power exponential correlation
             % type=2, cubic correlation, anisotropic
fittype=1;  % REML estimate

%%%% making F matrix
n = length(Y(:,1));
p = length(X(1,:));

%%%% Just try to use the loglik function.
if ctype==0,
    cor=zeros(2*p,1);
    cor(1:p) = corr;
    cor(p+1:2*p) = 2;
    loglik = llik1d(cor, Y, X, ctype, fittype, F, 2*p);
elseif ctype==1,
    loglik = llik1d(corr, Y, X, ctype, fittype, F, 2*p);
elseif ctype  ==  2,
    %% likelihood value for cubic correlation %%
    loglik = llik1d(corr, Y, X, ctype, fittype, F, 2*p);
end