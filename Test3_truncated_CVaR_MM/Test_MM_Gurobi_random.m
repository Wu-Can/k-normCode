%%=============================================================================================
%% Test the MM+Gurobi for the nonconvex truncated CVaR-based model with random data
%% Input:
%% m_vec = sample size vector;
%% n_vec = feature size vector;
%% lamc_vec = the vector of parameter lamc satisfying lambda:=lamc*(k_1-k_2).
%%============================================================================================

clear all; clc;
rng('default');
HOME = pwd;
addpath(genpath(HOME));
%%
%profile on
%================================== INPUT =================================
m_vec = 500;%[500, 500, 1000, 1500, 3000];
n_vec = 3000;%[3000, 50000, 100000, 200000, 500000];
lamc_vec = [0.15 0.1 0.05];
%==========================================================================
%% -------------------------Generate A and b---------------------
len_m = length(m_vec);
len_n = length(n_vec);
if len_m ~= len_n
   fprintf('m_vec and n_vec must be the same length !'); 
end
lenlam = length(lamc_vec);
result = zeros(len_m*lenlam,13);

for ii = 1:len_m
    m = m_vec(ii);
    n = n_vec(ii);
    [A,b] = Generate_A_b_Conf_large(m,n,2);
    
    %---------------------input parameter ---------------
    OPTIONS.tol = 1e-6;
    OPTIONS.m = m;
    OPTIONS.n = n;
    OPTIONS.rho0 = 1;
    OPTIONS.sigma0 = 1;
    OPTIONS.eta0 = 1;
    OPTIONS.rhoscale = 0.7;
    OPTIONS.sigmascale = 0.6;
    OPTIONS.etascale = 2;
    OPTIONS.rho_iter = 3;
    %--------------------------------------------------
    for jj = 1:lenlam
        kk1 = m;
        kk2 = 0.1*m;
        lamc = lamc_vec(jj);
        lambda = lamc*(kk1-kk2);
        
        [~,obj,info,runhist] = MM_Gurobi(A,b,kk1,kk2,lambda,OPTIONS);
        
        result(jj+(ii-1)*lenlam,1) = m;
        result(jj+(ii-1)*lenlam,2) = n;
        result(jj+(ii-1)*lenlam,3) = kk1;
        result(jj+(ii-1)*lenlam,4) = kk2;
        result(jj+(ii-1)*lenlam,5) = lamc;
        result(jj+(ii-1)*lenlam,6) = lambda;
        result(jj+(ii-1)*lenlam,7) = info.nnzeros_x;
        result(jj+(ii-1)*lenlam,8) = info.obj_gap;
        result(jj+(ii-1)*lenlam,9) = info.iter;
        result(jj+(ii-1)*lenlam,10) = info.iterBar;
        result(jj+(ii-1)*lenlam,11) = info.time;
        result(jj+(ii-1)*lenlam,12) = obj;
        
    end
end

save result_MM_Gurobi_random.mat result
%profile viewer

