%%=============================================================================================
%% Comparison of the truncated CVaR-based linear regression model and Average model
%% Input:
%% repeat_times = the number of repetitions of the simulation;
%% conf = 1, generate data by Configuration 1 under (m,n)=(100,1000);
%%        2, generate data by Configuration 2 under (m,n)=(100,20000);
%% lambdavec = the vector lambda;
%% num_kk2_vec = coefficients of k2;
%% flag_cm = 1, no contamination in sampling instance;
%%           2, vertical outliers in sampling instance;
%%           3, leverage points in sampling instance;
%%=============================================================================================
% profile on
clc;clear;
%================================= INPUT ==================================
repeat_times = 100;
conf = 1;
lambdavec = [1e-11 1e-10 1e-9 1e-8 1e-7 1e-6 1e-5 1e-4 1e-3 1e-2 1e-1 1:60]; 
num_kk2_vec = [0.15, 0.1, 0.05, 0];
flag_cm_vec = [1,2,3];
%==========================================================================
repeat_start = 1;
repeat_end = repeat_times;
lenkk2 = length(num_kk2_vec);
lenflag = length(flag_cm_vec);
result = zeros(lenkk2,8*lenflag);
Result = zeros(lenkk2,8*lenflag);
Time = zeros(repeat_times,1);
Oracle = 0;
tstart = clock;
for iii = repeat_start:repeat_end
    %% generate (A_valid,b_valid) and (A,b)
    rng(iii + repeat_times)
    [A_valid,b_valid,oracle] = Generate_A_b_Configuration(conf,1,1);
       for flag = 1:lenflag
            flag_cm = flag_cm_vec(flag);
            rng(iii)
            [A,b,~] = Generate_A_b_Configuration(conf,flag_cm,0);
            for prob = 1:lenkk2
                num_kk2 =  num_kk2_vec(prob);
                 result(prob,1+8*(flag-1):8*flag) = Cross_validation_truncated_CVaR...
                    (A,b,A_valid,b_valid,lambdavec,num_kk2,iii,flag_cm);
            end
       end
    Result = Result + result;
    Oracle = Oracle + oracle;
    ttime = etime(clock,tstart);
    Time(iii) = ttime;
end
Ave_Result = (1/repeat_times)*Result;
Ave_Oracle =  (1/repeat_times)*Oracle;
% 
save Ave_Result_nonconvex.mat Ave_Result
save Ave_Oracle_nonconvex.mat Ave_Oracle
%profile viewer


