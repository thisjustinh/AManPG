clear;
close all;
addpath ../../misc
addpath ../../SSN_subproblem
addpath SCCALab
addpath TFOCS-master

p = 300; % dim of X
q = 600; % dim of Y
N = 200; % sample size
r = 2;   % column number
type = 1; % matrix type

tic;  fid =1;
fprintf(fid,  'generating data.... \n');
% mean
mu_X = zeros(p,1);   mu_Y = zeros(q,1);
%% generate data matrix
% Identity
if type == 1
    sigma_X = eye(p);   sigma_Y = eye(q);
else
    if type == 2
        % toeplitz matrix
        a = 0.3;
        c1 = a.^((1:p)-1); c2 = a.^((1:q)-1);
        sigma_X = toeplitz(c1);   sigma_Y = toeplitz(c2);
    else
        % sparse inverse
        sigma_X = zeros(p); sigma_Y = zeros(q);
        c1 =zeros(p,1); c1(1:3,:) = [1; 0.5; 0.4]; omegaX = toeplitz(c1);    sigma_X0 = inv(omegaX);
        c2 =zeros(q,1); c2(1:3,:) = [1; 0.5; 0.4]; omegaY = toeplitz(c2);    sigma_Y0 = inv(omegaY);
        for i=1:p
            for j=1:p; sigma_X(i,j)= sigma_X0(i,j)/(sqrt(sigma_X0(i,i)*sigma_X0(j,j))); end
        end
        for i=1:q
            for j=1:q; sigma_Y(i,j)= sigma_Y0(i,j)/(sqrt(sigma_Y0(i,i)*sigma_Y0(j,j))); end
        end
        sigma_X(abs(sigma_X)<1e-3)=0;  sigma_Y(abs(sigma_Y)<1e-3)=0;
    end
end
% vector and rho
u1= zeros(p,r);  v1 = zeros(q,r);
s  = [1, 6, 11, 16, 21];
lambda = diag( [0.9;0.8] );
%rng('default')
%rng(10);
u1(s,(1:r)) = randi( [-2, 2], size(u1(s,(1:r))));

Tss = sigma_X(s, s);
u1 = u1 / sqrtm(u1(s,1:r)' * Tss * u1(s,1:r)+ 1e-13*eye(r));
%rng(10);
v1(s,(1:r)) = randi( [-2, 2], size(v1(s,(1:r))));
Tss = sigma_Y(s, s);
v1 = v1 / sqrtm(v1(s,1:r)' * Tss * v1(s,1:r)+ 1e-13*eye(r));


[u_n, ~, ~] = svd(u1, 0);
[v_n, ~, ~] = svd(v1, 0);
%generate covariance matrix
sigma_XY = sigma_X*u1*lambda*v1'*sigma_Y;

Data =  mvnrnd([mu_X; mu_Y] ,[sigma_X,sigma_XY; sigma_XY', sigma_Y],N); % data matrix
time = toc; fid =1;
fprintf(fid,  'generating data time uses %3.3f seconds\n',time );


%% run  scca_init to obtain initial point (you can use other intial points)
% convex relaxation 
% set parameter
b_set = [0.8;1;1.2;1.4;1.6];
% amanpg_set = [0.6;0.8;1;1.2];
manpg_set = 0.5*b_set;
time_manpg = 0;
%time_adap = 0;
time_init = 0;

Xtrain = Data(:,1:p);     Xtest = Xtrain;
Ytrain = Data(:,p+1:p+q); Ytest = Ytrain;
%[xhat,~] = scca_init(Xtrain, Ytrain, r,0.55, 1e-1, 1);% pre_run to initialize tfocs
%%init_stag;
tic;
[xhat,~] = scca_init(Xtrain, Ytrain, r,0.55, 1e-4, 1);
[U0, S0, V0] = svds(xhat, r);
Uinit_1 = U0(:,1:r);
Vinit_1 = V0(:,1:r);
[u0hat, ~,~] = svd(Uinit_1,'econ');  [v0hat,~,~] = svd(Vinit_1,'econ');
u0hat = u0hat(:,1:r);     v0hat = v0hat(:,1:r);
Init_lossu_1 = norm(u0hat * u0hat'  - u_n * u_n', 'fro')^2;
Init_lossv_1 = norm(v0hat * v0hat'  - v_n * v_n', 'fro')^2;
time_init_1 =toc;

%% AManPG
for cross_b = 1:5
    
    tic;
    manpg_para = manpg_set(cross_b);
    Uinit_manpg = Uinit_1;  
    Vinit_manpg = Vinit_1;

    option.b1 = manpg_para;%*sqrt((r+log(p))/N); 
    option.b2 = manpg_para;%*sqrt((r+log(q))/N); 
    option.maxiter = 1e3;    option.tol =1e-8;   option.inner_tol = 1e-10;
    option.n = r;  % column number
    option.q = q;
    option.p = p;
    option.X0 = Uinit_manpg;  %initial point
    option.Y0 = Vinit_manpg;  % initial point
    
    option.Xtest = Xtest;  option.Ytest = Ytest;
    
    %[xhat,~] = scca_init(Xtrain, Ytrain, r,option.tau1, 1e-4, 1);
    fprintf('=========================================================\n');

    fprintf('running A-Manpg............\n');
    %% amanpg_init_1

    [result_manpg_alt] = scca_amanpg(Xtrain,Ytrain,option,'l21');
    
    %% compute loss
    [uhat, ~,~] = svd(result_manpg_alt.X,0);  [vhat,~,~] = svd(result_manpg_alt.Y,0);
    result_manpg_alt.lossu = norm(uhat * uhat'  - u_n * u_n', 'fro')^2;
    result_manpg_alt.lossv = norm(vhat * vhat'  - v_n * v_n', 'fro')^2;
    % [~,~,Res.pho]  = canoncorr(A * uhat, B * vhat);
    if r == 1
        result_manpg_alt.lossu = 2*(1-abs(uhat'*option.u_n));
        result_manpg_alt.lossv = 2*(1-abs(vhat'*option.v_n));
    end
    [~,~,result_manpg_alt.rho]  = canoncorr(option.Xtest * uhat, option.Ytest * vhat);
    fprintf('------------------result of A-Manpg :--------------------\n');
   
    fprintf('Canonical correlations on test data:  rho = %2.3f; \n',result_manpg_alt.rho);
    fprintf('Projection U error: %2.3f; \n',result_manpg_alt.lossu);
    fprintf('Projection V error: %2.3f. \n',result_manpg_alt.lossv);
    fprintf('total_time: %.3f \n', result_manpg_alt.time);
    result_manpg_alt_1 = result_manpg_alt;
    time_manpg_1 = toc + time_init_1;
    

end