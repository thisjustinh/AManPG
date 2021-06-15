function [ Res] = scca_amanpg(A,B,option,sp_type)
%min -Tr(X'A'*B*Y)+ tau1*||X||_1 + tau2 ||Y||_1,  ------ sp_type = 'l1'
% or min -Tr(X'A'*B*Y)+ tau1*||X||_{2,1} +tau2 ||Y||_{2,1},  -- sp_type='l21'

%  s.t. X'*A'*A*X=I_n,  Y'*B'*B*Y=In
%% manpg alternating
%parameters

tic;
n=option.n;  % number of column
p=option.p;  %dim A
q=option.q;  %dim B
m = size(A,1);%sample
tau1 = option.b1*sqrt((n+log(p))/m);
tau2 = option.b2*sqrt((n+log(q))/m);
fprintf('A-ManPG parameter tau1: %2.3f; \n', tau1);
fprintf('A-ManPG parameter tau2: %2.3f; \n', tau2);

maxiter =option.maxiter;
tol = option.tol;
inner_tol =  option.inner_tol;
if strcmp(sp_type,'l1') 
    h1=@(X) tau1*sum(sum(abs(X)));
    h2=@(X) tau2*sum(sum(abs(X)));
    prox_func = @(b,lambda,r) proximal_l1(b,lambda,r);
end
if  strcmp(sp_type,'l21')
    h1=@(X) tau1*sum(vecnorm(X,2,2));
    h2=@(X) tau2*sum(vecnorm(X,2,2));
    prox_func = @(b,lambda,r) proximal_l21(b,lambda,r); 
end
inner_flag1=0;
%setduplicat_pduplicat(n);
Dn = sparse(DuplicationM(n));
pDn = (Dn'*Dn)\Dn';
%center data
A = A - repmat(mean(A,1),m,1);
B = B - repmat(mean(B,1),m,1);
t_min = 1e-4; % minimum stepsize
AtB = A'*B/(m-1); % covariance
BtA = AtB';
AtA = A'*A/(m-1);
BtB = B'*B/(m-1);
%L = max( svds(AtB,1), svds(AtB',1));  
%% set type
gamma = 1e-4/(m-1);
if svds(A,1,'smallest') < 1e-4
    M1 = (1-gamma)*AtA + gamma*eye(p); pdA = 0;
else
    M1 = AtA; pdA = 1;
end
if svds(B,1,'smallest') < 1e-4
    M2 = (1-gamma)*BtB + gamma*eye(q); pdB = 0;
else
    M2 = BtB;  pdB =1;
end
if m > p/2;    typeA = 1;   else;  typeA=0;  end  % p^2 n
if m > q/2;    typeB = 1;   else;  typeB = 0;  end
%% initial point
% X = option.X0;  Y= option.Y0;
X0 = option.X0; Y0 = option.Y0;
% [xhat,~] = scca_init(A, B, n, 0.55, 1e-4, 1);
% [U0, ~, V0] = svds(xhat, n);
X0 = X0(:,1:n);
Y0 = Y0(:,1:n);
%[U1,V1] = canoncorr(A,B); X0 = U1(:,1:n);  Y0 = V1(:,1:n);
X0 = retraction(X0,M1,A,typeA); 
Y0 = retraction(Y0,M2,B,typeB);
%time_init = toc;
X = X0;   Y = Y0;
% fprintf('------------------result of initial_point :--------------------\n');
%
% [uhat, ~,~] = svd(X,0);  [vhat,~,~] = svd(Y,0);
% Init_lossu = norm(uhat * uhat'  - option.u_n * option.u_n', 'fro')^2;
% Init_lossv = norm(vhat * vhat'  - option.v_n * option.v_n', 'fro')^2;
% [~,~,Init_rho]  = canoncorr(option.Xtest * uhat, option.Ytest * vhat);
%
% fprintf('Canonical correlations on test data:  rho = %2.3f;  \n',Init_rho);
% fprintf('Projection U error: %2.3f; \n',Init_lossu);
% fprintf('Projection V error: %2.3f. \n',Init_lossv);
% fprintf('time: %.3f \n', time_init);
%%
if m > p/2;   ABX = -(AtB'*X);  else;  ABX = -(B'*(A*X))/(m-1);   end  % p^2 n
if m > q/2;   ABY = -(AtB*Y);    else;  ABY = -(A'*(B*Y))/(m-1);  end

if typeA ==1;    MX = M1*X;
else
    if pdA ==0
        MX = (1-gamma)* A'*(A*X)/(m-1) + gamma*X;
    else ; MX = A'*(A*X)/(m-1);
    end
end
if typeB ==1;    MY = M2*Y;
else
    if pdB ==0
        MY = (1-gamma)* B'*(B*Y)/(m-1) + gamma*Y;
    else ; MY = B'*(B*Y)/(m-1);
    end
end

F(1) = sum(sum(Y.*(ABX))) + h1(X) + h2(Y);
num2 = zeros(maxiter,1); 
num1 = num2;
flag_maxiter = 0;% flag_linesearch =zeros(maxiter,1);
num_inex = 0; 
linesearch_num = 0;
alpha = 1;
normDsquared_Y =1;
fX = sum(sum(X.*(ABY))) + h1(X);
for iter=2:maxiter
    %% update X
    gx = ABY;  pgx=gx;  % grad or projected gradient both okay
    %% subproblem
    if alpha < t_min || num_inex >10
        inner_tol = option.inner_tol/100; % if subproblem inexact, decrease the tol
    else
        inner_tol = option.inner_tol;
    end
    t1 = 1;
    if iter == 2
        [ PX,num1(iter),Lam1,~,in_flag1]=Semi_newton_matrix_l21(p,n,MX,t1,X-t1*pgx,tau1*t1,inner_tol,prox_func,zeros(n),Dn,pDn);
    else
        [ PX,num1(iter),Lam1,~ ,in_flag1]=Semi_newton_matrix_l21(p,n,MX,t1,X-t1*pgx,tau1*t1,inner_tol,prox_func,Lam1,Dn,pDn);
    end
    if in_flag1 == 1   % subprolem total iteration.
        inner_flag1 = 1 + inner_flag1;
    end
    alpha=1;
    DX = PX - X; %descent direction D
    X_temp  = retraction(PX,M1,A,typeA);
    %fX = sum(sum(X.*(AY))) + h1(X);
    f_trial=sum(sum(X_temp.*(ABY)));
    f_trialX=f_trial+ h1(X_temp) ;   normDsquared_X=norm(DX,'fro')^2 ;
    if  max( normDsquared_X,  normDsquared_Y) < tol
        fprintf('A-ManPG terminates: converged iteration:%4d\n', iter);
        break;
    end
    %% linesearch
    while f_trialX >= fX - 1e-4*alpha*normDsquared_X
        alpha=0.5*alpha;
        if alpha < t_min
            num_inex = num_inex+1;
            break;
        end
        PX=X+alpha*DX;   X_temp  = retraction(PX,M1,A,typeA);
        linesearch_num = linesearch_num +1;
        f_trial = sum(sum(X_temp.*(ABY)));    f_trialX = f_trial + h1(X_temp) ;
    end
    
    X = X_temp;
    if m > p/2;   ABX = -(BtA*X);   else;  ABX = -(B'*(A*X)/(m-1));    end
    
    if typeA ==1;    MX = M1*X;
    else
        if pdA ==0
            MX = (1-gamma)* A'*(A*X)/(m-1) + gamma*X;
        else ; MX = A'*(A*X)/(m-1);
        end
    end
    %%  update Y
    gy = ABX;  pgy=gy;  t2 =1;
    if iter == 2
        [ PY,num2(iter),Lam2,~,in_flag2] = Semi_newton_matrix_l21(q,n,MY,t2,Y-t2*pgy,tau2*t2,inner_tol,prox_func,zeros(n),Dn,pDn);
    else
        [ PY,num2(iter),Lam2,~ ,in_flag2] = Semi_newton_matrix_l21(q,n,MY,t2,Y-t2*pgy,tau2*t2,inner_tol,prox_func,Lam2,Dn,pDn);
    end
    alpha=1;     DY = PY - Y;%descent direction D
    Y_temp  = retraction(PY,M2,B,typeB);
    
    fY = sum(sum(Y.*(ABX))) + h2(Y);
    f_trial=sum(sum(Y_temp.*(ABX)));
    f_trialY=f_trial+ h2(Y_temp) ;   normDsquared_Y = norm(DY,'fro')^2 ;
    if max( normDsquared_X,  normDsquared_Y) < tol
        fprintf('A-ManPG terminates: converged iteration:%4d\n', iter);
        
        break;
    end
    %% linesearch
    while f_trialY >= fY -1e-4*alpha*normDsquared_Y
        alpha=0.5*alpha;
        if alpha < t_min
            num_inex = num_inex+1;
            break;
        end
        PY = Y + alpha*DY;    Y_temp  = retraction(PY,M2,B,typeB);
        linesearch_num = linesearch_num +1;
        f_trial = sum(sum(Y_temp.*(ABX)));  f_trialY = f_trial+ h2(Y_temp) ;
    end
    Y = Y_temp;
    
    if m > q/2;   ABY = -(AtB*Y_temp);   else;  ABY = -(A'*(B*Y_temp)/(m-1));   end
    
    if typeB ==1;    MY = M2*Y;
    else
        if pdB ==0
            MY = (1-gamma)* B'*(B*Y)/(m-1) + gamma*Y;
        else ; MY = B'*(B*Y)/(m-1);
        end
    end
    fX = sum(sum(X.*(ABY))) + h1(X);
    F(iter) = fX + h2(Y);
    
    if iter ==maxiter
        flag_maxiter =1;
        fprintf('A-ManPG terminates: Achieved maximum iteration. \n');
    end
    
    
end
X((abs(X)<=1e-4))=0;
Y((abs(Y)<=1e-4))=0;
%X_manpg=X;  Y_manpg = Y;
%Res.init_lossu = Init_lossu;
%Res.init_lossv = Init_lossv;

Res.time = toc;
Res.sparsityX = sum(sum(X~=0));%/(p*n);
Res.sparsityY = sum(sum(Y~=0));%/(q*n);
Res.Fval =  F(iter-1);
Res.X = X; 
Res.Y =Y;
Res.flag_maxiter =flag_maxiter;  
Res.iter = iter;
Res.inner_total = sum(num1)+sum(num2); 
Res.linesearch_num = linesearch_num;

% [uhat, ~,~] = svd(X,0);  [vhat,~,~] = svd(Y,0);
% Res.lossu = norm(uhat * uhat'  - option.u_n * option.u_n', 'fro')^2;
% Res.lossv = norm(vhat * vhat'  - option.v_n * option.v_n', 'fro')^2;
% % [~,~,Res.pho]  = canoncorr(A * uhat, B * vhat);
% if n == 1
% Res.lossu = 2*(1-abs(uhat'*option.u_n));
% Res.lossv = 2*(1-abs(vhat'*option.v_n));
% end
% [~,~,Res.rho]  = canoncorr(option.Xtest * uhat, option.Ytest * vhat);
%clear F; %clear num2;
%clear global Dn pDn
end