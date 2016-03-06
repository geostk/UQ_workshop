function cor = correxp(Xpred,X,theta,pow);
%%%% Calculate the product power exponential correlation matrix, where
%%%% for i=1:npred and j=1:n,
%%%% cor(i,j) = exp{-sum_(from k=1 to k=p){theta(k)*(xpred(i,k)-x(j,k))^pow(k)}}.
%%%% given Xpred: npred by p; X: n by p; theta: p by 1, and pow: p by 1, 

[n,p] = size(X); %%%% check the size %%%%
[npred,p] = size(Xpred);
cor = ones(npred,n);
p1 = length(theta);

if p1 ~= p,
    fprintf('Dimensions of X and theta are not the same');
else
%%%% Making the covariance matrix %%%%
    for i = 1:npred,
        for j = 1:n,
            for k = 1:p,
                cor(i,j) = cor(i,j) * exp(-theta(k)*abs(Xpred(i,k)-...
                    X(j,k))^pow(k)); 
            end
        end
    end
end     
