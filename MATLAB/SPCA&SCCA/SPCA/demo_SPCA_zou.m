%function compare_spca

clear% all;
close all;
addpath ../../misc
addpath ../../SSN_subproblem
d = 500; % dimension
m = 1000; %sample size
n = 4; % column number
mupar=[0.1;0.2;0.3;0.4;0.5;0.6]; % sparsity parameters

fid =1;
for k = 1:2 %mu
    lambda = 1;
    mu = mupar(k);
    fprintf(fid,'========================================================\n');
    
    for test_random =1:1   %multiple times ----average. 
        fprintf(fid,' d:%4d, m:%4d, n:%3d, mu:%2.2f, lambda:%2.2f \n',d,m,n,mu,lambda);

        %fprintf(fid,'%4d %4d %3d %2.2f \n',d,m,n,lambda);
        fprintf(fid,'--iter----Fval----sparsity----cpu\n');
        %%%%%% random testing
        rng('default');
        rng(10);
        A = randn(m,d);     mu = mupar(k)*ones(n,1)';  type =0;
        if (type == 1) %covariance matrix
            scale = max(diag(A)); % Sigma=A/scale;
        elseif (type == 0) %data matrix
            A = A - repmat(mean(A,1),m,1);
            scale = [];
            for i = 1:d
                scale = [scale norm(A(:,i))];
            end
            scale = max(scale);
            A =A/scale;
            A1 = A;   AtA = A'*A;
        end
        %%%%%% random intialization
        %    rng(10);
        % [phi_init,~] = svd(randn(d,n),0);
        
        %%%%%% singular vectors initialization
        [~,~,Vu] = svds(A,n);  
        phi_init = Vu;   
        
        %%%%%% AManPG
        option_manpg_alt.X0 = phi_init;
        option_manpg_alt.maxiter =1e4; 
        option_manpg_alt.tol =1e-5;
        option_manpg_alt.n = n;  
        option_manpg_alt.d = d;  
        option_manpg_alt.mu = mu;  
        option_manpg_alt.type = type;
        option_manpg_alt.lambda = lambda;
        
        [iter_manpg_alt(test_random), F_man_alt(test_random),sparsity_alt(test_random),time_alt(test_random), Y_alt]= spca_amanpg(A,option_manpg_alt);
        
        print_format =  '%5d  %1.3e   %1.3f    %3.2f \n';
        fprintf(fid,print_format, mean(iter_manpg_alt),    mean( F_man_alt),mean(sparsity_alt),mean(time_alt));
 
    end
end





