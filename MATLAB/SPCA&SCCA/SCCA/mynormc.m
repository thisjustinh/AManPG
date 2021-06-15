function X = mynormc(X)

n = size(X,2); 
for j=1:n
    normj = norm(X(:,j));
    X(:,j) = X(:,j)/normj;
end