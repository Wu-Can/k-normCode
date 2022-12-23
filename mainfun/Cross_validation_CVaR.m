%%*************************************************************************
%% Cross_validation_CVaR: five-fold cross validation according to 
%% TMAPE(lambda) to select the parameter lambda of the following 
%% convex CVaR-based sparse linear regression model
%%
%% (P)   minimize_{x in R^n} {||A*x - b||_(k) + lambda ||x||_1}
%%
%% solution_vec = Cross_validation_CVaR(A,b,A_valid,b_valid,lambdavec,kk1,iii)
%%
%% Input: 
%% A, b = matrix A and vector b in (P)
%% A_valid, b_valid = an additional validation set
%% lambdavec = the vector of the values of parameter lambda in (P)
%% kk1 = the value of parameter k in (P)
%% iii = the number of repeat times
%% Output:
%% solution_vec = result corresponding to the minimal value of 
%%                TMAPE(lambda) defined in Appendix F.1 of the 
%%                supplementary materials
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the Appendix F.1 of the 
%% supplementary materials in the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%*************************************************************************
function solution_vec = Cross_validation_CVaR(A,b,A_valid,b_valid,lambdavec,kk1,iii)
m_valid = length(b_valid);
num_df = 5;
lenlam = length(lambdavec);
result = zeros(lenlam*num_df,8);
Result = zeros(lenlam,8);

[m,n] = size(A);
m_part = m/num_df;
index_true_nnz = [1:11];
index_true_0 = [12:n];
len_index_true_nnz = length(index_true_nnz);
len_index_true_0 = length(index_true_0);

for jj = 1:lenlam
    OPTIONS.lambda = lambdavec(jj);
    for ii = 1:num_df
        fprintf('*-----------------------------------------* \n');
        fprintf('Repeat=%3.0d, kk1=%3.2f, lambda=3.2e%, Cross_validation_DC=%3.0d',iii,kk1,OPTIONS.lambda,ii);
        A_train = A;
        b_train = b;
        A_test = A(1+m_part*(ii-1):ii*m_part,:);
        b_test = b(1+m_part*(ii-1):ii*m_part);
        A_train(1+m_part*(ii-1):ii*m_part,:) = [];
        b_train(1+m_part*(ii-1):ii*m_part) = [];
        
        m_test = length(b_test);
        [m_train,n_train] = size(A_train);
        OPTIONS.tol = 1e-6;
        OPTIONS.flag_tol = 0;
        OPTIONS.m = m_train;
        OPTIONS.n = n_train;
        
        OPTIONS.sigma = 20/OPTIONS.lambda;
        OPTIONS.sigmascale = 1.25;
        OPTIONS.tau = 20/OPTIONS.lambda;
        OPTIONS.tauscale = 1.2;
        if kk1 == 0
            OPTIONS.kk = 1;
            kk_test = 1;
        else
            OPTIONS.kk = ceil(kk1*m_train);
            kk_test = ceil(kk1*m_test);
        end
        
        [~,~,~,~,~,info] = NALM(A_train,b_train,OPTIONS); 
        Ax_testb = A_test*info.x-b_test;
        Ax_testb_tmp = sort(abs(Ax_testb),'descend');
        Ax_valid = A_valid*info.x;
        MAPE = norm(Ax_valid-b_valid,1)/m_valid;
        TMAPE = (1/kk_test)*sum(Ax_testb_tmp(1:kk_test));
        FPR = length(intersect(info.index_nnz,index_true_0))/len_index_true_0;
        FNR = length(intersect(info.index_0,index_true_nnz))/len_index_true_nnz;

        result(ii+num_df*(jj-1),1) = OPTIONS.lambda;
        result(ii+num_df*(jj-1),2) = info.xnnz;
        result(ii+num_df*(jj-1),3) = MAPE;
        result(ii+num_df*(jj-1),4) = FPR;
        result(ii+num_df*(jj-1),5) = FNR;        
        result(ii+num_df*(jj-1),6) = norm(info.x,1);        
        result(ii+num_df*(jj-1),7) = info.eta_res;
        result(ii+num_df*(jj-1),8) = TMAPE;
         
        if (OPTIONS.lambda > 1e-6) && (ii == 1)
            OPTIONS.x0 = info.x;
            OPTIONS.z0 = info.z;
            OPTIONS.u0 = info.u;
        end
    end
    Result(jj,:) = mean(result(1+num_df*(jj-1):num_df+num_df*(jj-1),:));
end

[~, minindex] = min(Result(:,end));
solution_vec = Result(minindex,:);

