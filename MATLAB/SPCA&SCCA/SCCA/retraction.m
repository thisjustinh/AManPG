function [X ]= retraction(Y,M,B,type)
% retraction  % modified from manopt toolbox  P.-A Absil.
% Y : tangent point;            M: generalized Stiefel manifold
% B : M = BTB;
% type = 1;   using covariance matrix
% otherwise   data matrix preferred
m = size(B,1);
[u, ~, v] = svd(Y, 0);
% Instead of the following three steps, an equivalent, but an
% expensive way is to do X = u*(sqrtm(u'*(B*u))\(v')).
if type == 1
    [q, ssquare] = eig(u'*(M*u));
else
    [q, ssquare] = eig(u'*(B'*(B*u))/(m-1));
end
qsinv = q/sparse(diag(sqrt(diag(ssquare))));
X = u*((qsinv*q')*v'); % X'*B*X is identity.
