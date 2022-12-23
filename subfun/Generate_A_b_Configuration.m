%%=========================================================================
%% Generate_A_b_Configuration:
%% generate the small-scale random data A and b by the following 
%% linear regression model:
%%                     b = A*x_true + epsilon,
%% where A is a design matrix in R^{m*n}, b is the response vector
%% in R^m, epsilon is an error vetor in R^m.
%%
%% [A,b,oracle] = Generate_A_b_Configuration(conf,flag_cm,orc)
%%
%% Input: 
%% conf = 1, generate random data by Configuration 1 (m=100, n=1000);
%%        2, generate random data by Configuration 2 (m=100, n=20000);
%% flag_cm = 1, no contamination;
%%           2, vertical outliers;
%%           3, leverage points;
%% orc = 1, calculate oracle:=(1/m)||A*x_true - b||_1; 
%%          otherwise, does not calculate 
%% Output: A, b, oracle
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the Table 3 in the supplementary 
%% materials of the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%=========================================================================
function [A,b,oracle] = Generate_A_b_Configuration(conf,flag_cm,orc)
m = 100;
outnum = 0.1*m; % the number of outliers
oracle = 0;
switch(conf)
    case 1
        n = 1000; a = 0.5;
        %% Generate A
        A = mvnrnd(zeros(1,n),a.^toeplitz(0:n-1),m);
        if (flag_cm == 3)
            A_tmp = 50 + randn(outnum,n);
            outpos = randperm(m,outnum);
            A(outpos,:) = A_tmp;
        end
        %% Generate x_true
        x_true = zeros(n,1);
        x_true([1,7]) = 1.5;
        x_true(2) = 0.5;
        x_true([4,11]) = 1;
        Ax_ture = A*x_true;
        %% Generate epsilon and b
        epsilon = 0.5*randn(m,1);
        if (flag_cm == 2) || (flag_cm == 3)
            epsilon_tmp = 20 + 0.5*randn(outnum,1);
            outpos = randperm(m,outnum);
            epsilon(outpos) = epsilon_tmp;
        end
        b = Ax_ture + epsilon;
    case 2
        n = 20000; a = 0.6;
        %% Generate matrix A
        n1 = 1000; n2 = n - n1;
        A_tmp1 = mvnrnd(zeros(1,n1),a.^toeplitz(0:n1-1),m);
        A_tmp2 = randn(m,n2);
        A = [A_tmp1,A_tmp2];
        if (flag_cm == 3)
            A_tmp = 50 + randn(outnum,n);
            outpos = randperm(m,outnum);
            A(outpos,:) = A_tmp;
        end
        %% Generate x_true
        x_true = zeros(n,1);
        x_true([1:10]) = 1;
        Ax_ture = A*x_true;
        %% Generate epsilon and b
        epsilon = randn(m,1);
        if (flag_cm == 2) || (flag_cm == 3)
            epsilon_tmp = 20 + randn(outnum,1);
            outpos = randperm(m,outnum);
            epsilon(outpos) = epsilon_tmp;
        end
        b = Ax_ture + epsilon;
end
if orc == 1
    oracle = (1/m)*norm(epsilon,1);
end










