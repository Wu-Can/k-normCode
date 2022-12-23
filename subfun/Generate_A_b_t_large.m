%%=========================================================================
%% Generate_A_b_t_large:
%% generate the large-scale random data A and b by the following 
%% linear regression model:
%%                      b = A*x_true + epsilon,
%% where A is a design matrix in R^{m*n}, b is the response vector 
%% in R^m,the error epsilon follows t(5), t(3) and Laplace(0,1), 
%% respectively. 
%%
%% [A,b] = Generate_A_b_t_large(m,n,flag_err)
%%
%% Input: 
%% m = the sample size
%% n = the feature size
%% flag_err = 1, epsilon follows trnd(5,m,1);
%%            2, epsilon follows trnd(3,m,1);
%%            3, epsilon follows eps = Laplace(0,1).
%% Output: A, b
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the Table 1 in the supplementary 
%% materials of the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%=========================================================================
function [A,b] = Generate_A_b_t_large(m,n,flag_err)

rng default;
a = 0.5;
%% generate A
n1 = n - 1;
T = Generate_Toeplitz_matrix(a,n1);
[R,idef] = chol(T);
if (idef)
    error('T is not positive definite');
end
A_tmp = randn(m,n1)*R;
A = [ones(m,1),A_tmp];

%% generate x_true
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
if flag_err == 1
    free = 5;
    epsilon = trnd(free,m,1);
elseif flag_err == 2
    free = 3;
    epsilon = trnd(free,m,1);
elseif flag_err == 3
    epsilon = randl(m,1);
end
b = Ax_ture + epsilon;

% eval(['save A_b_',num2str(m),'_',num2str(n),'_t',num2str(free),'.mat A b']);










