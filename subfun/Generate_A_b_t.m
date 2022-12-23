%%=========================================================================
%% Generate_A_b_t:
%% generate the small-scale random data A and b by the following 
%% linear regression model:
%%                             b = A*x_true + epsilon,
%% where A is a design matrix in R^{m*n}, b is the response vector 
%% in R^m,the error epsilon follows t(5), t(3) and Laplace(0,1), 
%% respectively. 
%%
%% [A,b,oracle] = Generate_A_b_t(m,n,flag_err,test)
%%
%% Input: 
%% m = the sample size
%% n = the feature size
%% flag_err = 1, epsilon follows trnd(5,m,1);
%%            2, epsilon follows trnd(3,m,1);
%%            3, epsilon follows eps = Laplace(0,1).
%% test = 1, calculate oracle:=(1/m)||A*x_true - b||_1; 
%%           otherwise, does not calculate 
%% Output: A, b, oracle
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the Table 1 in the supplementary
%% materials of the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%=========================================================================
function [A,b,oracle] = Generate_A_b_t(m,n,flag_err,test)
a = 0.5;
oracle = 0;
%% generate A
n1 = n - 1;
A_tmp = mvnrnd(zeros(1,n1),a.^toeplitz(0:n1-1),m);
A = [ones(m,1),A_tmp];

%% generate ground truth x_true
x_true = zeros(n,1);
for iii = 1:11
    if iii == 1
        x_true(iii) = 1;
    else
        x_true(iii) = 1/(iii-1);
    end
end

%% generate epsilon and b
Ax_ture = A*x_true;
switch(flag_err)
    case 1
        epsilon = trnd(5,m,1);
    case 2
        epsilon = trnd(3,m,1);
    case 3
        epsilon = randl(m,1);
end
b = Ax_ture + epsilon;

if test == 1
    oracle = (1/m)*norm(epsilon,1);
end






