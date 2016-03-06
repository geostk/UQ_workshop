function cor = cormatexp(X,theta,pow);
%%%% Calculate the product power exponential correlation matrix, where
%%%% cor(i,j) = exp{-sum_(from k=1 to k=p){theta(k)*(x(i,k)-x(j,k))^pow(k)}}.
%%%% given X == n by p; theta == p by 1, and pow == p by 1

[n p] = size(X); %%%% check the size %%%%
cor = ones(n,n);
p1 = length(theta);

if p1 ~= p
    fprintf('Dimensions of X and theta are not the same aa');
else
%%%% Making the covariance matrix.
    for i = 1:n,
        for j = 1:n,
            for k = 1:p,
                cor(i,j) = cor(i,j) * exp(-theta(k)*abs(X(i,k)-X(j,k))^pow(k)); 
            end
        end
    end
end     
