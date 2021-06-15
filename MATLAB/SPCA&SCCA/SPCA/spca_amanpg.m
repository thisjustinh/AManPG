function [iter,F_amanpg, sparsity, time,X,Y_man]= spca_amanpg(B,option,F_palm)
% Alternating rgd and pgd for spca, i.e.
% if option.lambda < inf,  min Tr(Y'*B'B*Y)- 2Tr(X'*B'*B*Y) + lambda*norm(Y,'fro')^2  + mu*norm(Y,1) s.t. X'*X=In.
% if option.lambda = inf,  min - 2Tr(X'*B'*B*Y)+ norm(Y,'fro')^2  + mu*norm(Y,1) s.t. X'*X=In.
% If option.type =0; B is a m*d matrix.....m is sample size.
% d is the number of feature.  B should be centered and normalized such that
% the column means are 0 and the column Euclidean lengths are 1.
% Or option.type = 1, B is Gram matrix, i.e. covariance matrix ;
% data matrix is suppose to be normalized.
% e.g.
% if (type == 1) %covariance matrix
%     scale = max(diag(A)); % Sigma=A/scale;
% elseif (type == 0) %data matrix
%     A = A - repmat(mean(A,1),m,1);
%     scale = [];
%     for i = 1:d
%         scale = [scale norm(A(:,i))];
%     end
%     scale = max(scale);
%     A =A/scale;
%     %  Sigma=A'*A;
% end
tic;
%parameters
mu = option.mu; 
lambda = option.lambda; 
n = option.n; 
type = option.type;
[m,d] = size(B);
%h=@(X) mu*sum(sum(abs(X)));
h=@(X) sum(mu.*sum(abs(X)));
maxiter =option.maxiter;
tol = option.tol;
if d < m*2
    B = B'*B;
    type = 1;
end
if type == 0
    LY = 2*(svds(B,1))^2 + 2*lambda;
else
    LY = 2*svds(B,1) + 2*lambda;
end

%% initial point
X0 = option.X0;  Y0= option.Y0;
X = X0;  Y = Y0;
total_linesearch = 0;    linesearch_flag = 1;  min_step = 0;
linesearch_flag_Y = 1;
t = 1/LY;     tau = 100/d;
if type == 0
    AY = B'*(B*Y);    AX = B'*(B*X);
else
    AY = B*Y;     AX = B*X;
end

if lambda ~= inf %elastic net parameter isn't inf
    fx =  -2*sum(sum(X.*(AY))); 
    fy = sum(sum(Y.*(AY)))+lambda*norm(Y,'fro')^2  + h(Y);
    F_rgd(1)= fx + fy;
    for iter = 2:maxiter
 
        % update Y
        if linesearch_flag_Y == 0
            t = t*1.01;
        else
            t = max(1/LY, t/1.01);
        end
         linesearch_flag_Y = 0;
        Y_t = proximal_l1(Y-t*2*(AY-AX+lambda*Y), mu*t,n);
        if type == 0; AYt = B'*(B*Y_t);   else;  AYt = B*Y_t;     end
        f_ytrial = -2*sum(sum(X.*(AYt))) + sum(sum(Y_t.*AYt)) + lambda*norm(Y_t,'fro')^2 + h(Y_t);
        normpg = norm(Y_t - Y,'fro')^2/t^2;
        while f_ytrial > F_rgd(iter-1) - 1e-3*t* normpg
            t = 0.5*t;
            if t < 1e-5/d
                break;
            end
            Y_t = proximal_l1(Y-t*2*(AY-AX+lambda*Y), mu*t,n);
             if type == 0; AYt = B'*(B*Y_t);   else;  AYt = B*Y_t;     end
            f_ytrial = -2*sum(sum(X.*(AYt))) + sum(sum(Y_t.*AYt)) + lambda*norm(Y_t,'fro')^2 + h(Y_t);
            linesearch_flag_Y = 1;
           
        end
       Y = Y_t;  
       AY = AYt;
        % update X
        if linesearch_flag == 0
            tau = tau*1.1;
        end
        if min_step == 1
            tau = 1/d;
        end
        linesearch_flag = 0;
        min_step = 0;
        gx = -2*AY;   xgx = gx'*X;    
        %RGX = gx - X*xgx;  %canonical riemannian gradient
        RGX = gx-0.5*X*(xgx+xgx');%projected gradient
        TX = X - tau*RGX;  
        [U, SIGMA] = eig(TX'*TX);  
        SIGMA = diag(SIGMA);  J =U*diag(sqrt(1./SIGMA))*U';
        X_trial= TX*J;    
        f_xtrial = - 2*sum(sum(X_trial.*AY));
        fXval = - 2*sum(sum(X.*AY));
        normpg = norm(RGX,'fro')^2;
        while f_xtrial > fXval - 1e-3*tau* normpg
            tau = 0.5*tau;
            if tau < 1e-5/d
                min_step = 1;
                break;
            end
            TX = X - tau*RGX;
            %  [U,~,V]=svd(TX, 0); X_trial=U*V';
            %  [X_trial,R]=qr(TX,0);      X_trial =X_trial*diag(sign(diag(R)));
            [U, SIGMA] = eig(TX'*TX);  
            %[U, SIGMA] = eig(eye(n)+ tau^2*GtG);  
            SIGMA = diag(SIGMA);  J =U*diag(sqrt(1./SIGMA))*U';
            X_trial= TX*J;
            total_linesearch = total_linesearch+1;
            linesearch_flag = 1;
            f_xtrial = - 2*sum(sum(X_trial.*AY));
        end
        
        X = X_trial;
        %  AX = B'*(B*X);
        if type == 0; AX = B'*(B*X);   else;  AX = B*X;     end
        
        fx = f_xtrial;%- 2*sum(sum(X.*AY));
        fy = sum(sum(Y.*AY)) + lambda*norm(Y,'fro')^2 + h(Y);
        F_rgd(iter)= fx + fy;
        
        %  if  normDsquared<tol^2
        if iter >1
            if  ( abs(F_rgd(iter)-F_rgd(iter-1))<tol && F_rgd(iter) < F_palm) || abs(F_rgd(iter)-F_rgd(iter-1))< 1e-12
                break;
            end
        end
    end
else%elastic net parameter is inf
    fx =  -2*sum(sum(X.*(AY))); fy = norm(Y,'fro')^2  + h(Y);
    F_rgd(1)= fx + fy;
    for iter = 2:maxiter
        if linesearch_flag == 0
            tau = tau*1.1;
        end
        if min_step == 1
            tau = 1/d;
        end
        min_step = 0;
        linesearch_flag = 0;
        % update Y
        %Y_old = Y;
        t = 1/2;
        Y = proximal_l1(Y-2*t*(-AX+ Y), mu*t,n);
        if type == 0;  AY = B'*(B*Y);   else;  AY = B*Y;     end
        % update X
        gx = -2*AY;   xgx = gx'*X;     RGX = gx - X*xgx;  %canonical riemannian gradient
        TX = X - tau*RGX;
        [U,~,V]=svd(TX, 0);   X_trial=U*V';
        f_xtrial = - 2*sum(sum(X_trial.*AY));
        fXval = - 2*sum(sum(X.*AY));
        normpg = norm(RGX,'fro')^2;
        while f_xtrial > fXval - 1e-3*tau* normpg
            tau = 0.5*tau;
            if tau < 1e-3/d
                min_step = 1;
                break;
            end
            TX = X - tau*RGX;
            [U,~,V]=svd(TX, 0); X_trial=U*V';
            total_linesearch = total_linesearch+1;
            linesearch_flag = 1;
            f_xtrial = - 2*sum(sum(X_trial.*AY));
        end
        
        X = X_trial;
        %  [U,~,V]=svd(2*AY,0); X=U*V'; minimization
        %  AX = B'*(B*X);
        if type == 0; AX = B'*(B*X);   else;  AX = B*X;     end
        
        fx = f_xtrial;%- 2*sum(sum(X.*AY));
        fy =  norm(Y,'fro')^2 +h(Y);
        F_rgd(iter)= fx + fy;
        
        %  if  normDsquared<tol^2
        if iter >1
            if  abs(F_rgd(iter)-F_rgd(iter-1))<tol
                %  error(iter) =normpg +(norm(Y -Y_old,'fro')^2/t2^2);
                % if normpg + norm(Y -Y_old,'fro')^2/t^2< tol
                break;
            end
        end
    end
    
end
F_amanpg = F_rgd(iter);
Y_norm = sqrt(sum(Y.^2));
Y_norm(Y_norm == 0) = 1;
Y_man = Y./(ones(d,1)*Y_norm);
%X_man = X;
%Y_man = Y;
time = toc;
sparsity= sum(sum(Y==0))/(d*n);
%fprintf('%1.5f\n', tau);
%semilogy(1:iter, F_rgd - F_manpg);
