%%=========================================================================
%% Comparison of the convex CVaR-based linear regression model,
%% Average model and Maximum model
%% Input: 
%% m = sample size;
%% n = feature size;
%% flag_err = 1, epsilon follows trnd(5,m,1);
%%            2, epsilon follows trnd(3,m,1);
%%            3, epsilon follows eps = Laplace(0,1).
%% repeat_times = the number of repetitions of the simulation;
%% num_kk_vec = the vecor of coefficients of k;
%% lambdavec = the vector lambda;
%%=========================================================================
%profile on
clc;clear;
%=============================== INPUT ====================================
m = 100;
n = 1000;
flag_err = 1; 
repeat_times = 100;
num_kk_vec = [1:-0.1:0];
lambdavec = [1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1:40];
%==========================================================================
repeat_start = 1;
repeat_end = repeat_times;
lenkk = length(num_kk_vec);
result = zeros(lenkk,8);
Result = zeros(lenkk,8);
Time = zeros(repeat_times,1);
Oracle = 0;
tstart = clock;
for iii = repeat_start:repeat_end
    %% Generate (A,b) and (A_valid,b_valid) 
    rng(iii + repeat_times)
    [A_valid,b_valid,oracle] = Generate_A_b_t(m,n,flag_err,1);
    rng(iii)
    [A,b,~] = Generate_A_b_t(m,n,flag_err,0);
    for prob = 1:lenkk
        num_kk =  num_kk_vec(prob);
        result(prob,:) = Cross_validation_CVaR(A,b,A_valid,b_valid,lambdavec,num_kk,iii);
    end
    Result = Result + result;
    Oracle = Oracle + oracle;
    ttime = etime(clock,tstart);
    Time(iii) = ttime;
end
Ave_Result = (1/repeat_times)*Result;
Ave_Oracle =  (1/repeat_times)*Oracle;

save results
save Ave_Result_convex.mat Ave_Result
save Ave_Oracle_convex.mat Ave_Oracle
%profile viewer




