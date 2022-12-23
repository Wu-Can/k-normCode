%%==========================================================================
%% Generate_Toeplitz_matrix:
%% generate n-dimensional Toeplitz covariance matrix T with
%%                      T(i,j)=a^(|i-j|)
%%
%% T = Generate_Toeplitz_matrix(a,n)
%%
%% Input:
%% a = the parameter a in the matrix T
%% n = the dimension of the matrix T
%% Output: T
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the supplementary materials of the
%% paper: Convex and Nonconvex Risk-based Linear Regression at Scale.
%%==========================================================================
function T = Generate_Toeplitz_matrix(a,n)
a_mul = 1;
for j = 1:n
    a_mul = a_mul*a;
    if a_mul < eps
        nnz_x = j;
        break;
    end
    nnz_x = j;
end
nnz_x = min(nnz_x,n-1);

x_tmp = zeros(nnz_x,1);
for i =1:nnz_x
    x_tmp(i) = a^i;
end

x = zeros(n,1);
if n <= nnz_x
    x = x_tmp(1:n);
else
    x(1:nnz_x) = x_tmp;
end

T_tmp = sparse(n,n);

for jj = 1:n-1
    t = min(n,nnz_x+jj);
    T_tmp(jj+1:t,jj) = x(1:t-jj);
end
T = T_tmp + speye(n) + T_tmp';
end

