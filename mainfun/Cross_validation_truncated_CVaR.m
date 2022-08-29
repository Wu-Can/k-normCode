%%
%% five-fold cross validation for solving the truncated CVaR-based problems on random data
%%
function solution_vec = Cross_validation_truncated_CVaR(A,b,A_valid,b_valid,lambdavec,num_kk2,iii,flag_cm)
m_valid = length(b_valid);
num_df = 5;

lenlam = length(lambdavec);
result = zeros(lenlam*num_df,8);
Result = zeros(lenlam,8);

[m,n] = size(A);
if n == 1000
    index_true_nnz = [1 2 4 7 11];
    index_true_0 = setdiff([1:n],index_true_nnz);
    len_index_true_nnz = length(index_true_nnz);
    len_index_true_0 = length(index_true_0);
elseif n == 20000
    index_true_nnz = [1:10];
    index_true_0 = [11:20000];
    len_index_true_nnz = length(index_true_nnz);
    len_index_true_0 = length(index_true_0);
end

m_part = m/num_df;
for jj = 1:lenlam
    lambda = lambdavec(jj);
    for ii = 1:num_df
        fprintf('Repeat: %3.0d, flag_cm: %3.0d, Cross_validation_DC: %3.0d',iii,flag_cm,ii);
        A_train = A;
        b_train = b;
        A_test = A(1+m_part*(ii-1):ii*m_part,:);
        b_test = b(1+m_part*(ii-1):ii*m_part);
        A_train(1+m_part*(ii-1):ii*m_part,:) = [];
        b_train(1+m_part*(ii-1):ii*m_part) = [];
        
        m_test = length(b_test);
        [m_train,n_train] = size(A_train);
        OPTIONS.tol = 1e-6;
        OPTIONS.m = m_train;
        OPTIONS.n = n_train;
        if num_kk2 > 0
            %% MM_NPPA for solving the nonconvex truncated CVaR-based model
            OPTIONS.Warm_starting = 0;
            OPTIONS.UCI = 0;
            OPTIONS.rho0 = 1;
            OPTIONS.sigma0 = 1;
            OPTIONS.eta0 = 1;
            OPTIONS.rhoscale = 0.7;
            OPTIONS.sigmascale = 0.6;
            OPTIONS.etascale = 2;
            OPTIONS.rho_iter = 3;
            kk1 = m_train;
            kk2 = ceil(num_kk2*m_train);
            kk1_test = m_test;
            kk2_test = ceil(num_kk2*m_test);
            
            [~,~,info,~] = MM_NPPA(A_train,b_train,kk1,kk2,lambda,OPTIONS);
            
            Ax_testb = A_test*info.x-b_test;
            Ax_testb_tmp = sort(abs(Ax_testb),'descend');
            TMAPE = (1/(kk1_test-kk2_test))*(sum(Ax_testb_tmp(1:kk1_test))-sum(Ax_testb_tmp(1:kk2_test)));
        else
            %% NALM for solving the convex CVaR-based model
            OPTIONS.lambda = lambda;
            OPTIONS.flag_tol = 0;
            OPTIONS.sigma = 20/OPTIONS.lambda;
            OPTIONS.sigmascale = 1.25;
            OPTIONS.tau = 20/OPTIONS.lambda;
            OPTIONS.tauscale = 1.2;
            OPTIONS.kk = m_train;
            
            [~,~,~,~,~,info] = NALM(A_train,b_train,OPTIONS);
            
            Ax_testb = A_test*info.x-b_test;
            TMAPE = norm(Ax_testb,1)/m_test;
        end
        Ax_valid = A_valid*info.x;
        MAPE = norm(Ax_valid-b_valid,1)/m_valid;
        FPR = length(intersect(info.index_nnz,index_true_0))/len_index_true_0;
        FNR = length(intersect(info.index_0,index_true_nnz))/len_index_true_nnz;

        result(ii+num_df*(jj-1),1) = MAPE;
        result(ii+num_df*(jj-1),2) = FPR;
        result(ii+num_df*(jj-1),3) = FNR;
        if num_kk2 > 0
            result(ii+num_df*(jj-1),4) = info.nnzeros_x;
            result(ii+num_df*(jj-1),5) = info.obj_gap;
        else
            result(ii+num_df*(jj-1),4) = info.xnnz;
            result(ii+num_df*(jj-1),5) = info.eta_res;
        end
        result(ii+num_df*(jj-1),6) = norm(info.x,1);
        result(ii+num_df*(jj-1),7) = lambda;
        result(ii+num_df*(jj-1),8) = TMAPE;      
    end
    Result(jj,:) = mean(result(1+num_df*(jj-1):num_df+num_df*(jj-1),:));
end
[~, minindex] = min(Result(:,end));
solution_vec = Result(minindex,:);


