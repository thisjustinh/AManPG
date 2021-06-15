function [vars_adj_rate] = vars_adj(u, A,AtA)
%adjusted variance ---- Zou 2003
 B = u'*AtA*u;% trace(u'*A*u); 
% [U D V] = svd(B);
 %Z1 = U*D.^0.5; 
Z = A*u;
[Q, R] = qr(Z);
%[~,R1] = qr(Z1);
r = diag(R);
var_adj = r.^2;
vars_adj_rate = var_adj./trace(AtA);