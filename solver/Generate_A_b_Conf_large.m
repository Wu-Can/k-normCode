%% Generate A and b by the following linear regression model:
%%                  b = A*x_true + epsilon,
%% where A is a design matrix in R^{m*n}, b is the response vector in R^m,
%% epsilon is a error vecor in R^m.
%% [A,b] = Generate_A_b_Conf_large(m,n,flag_cm)
%% Input: 
%% m = the sample size;
%% n = the feature size;
%% flag_cm = 1, no contamination;
%%           2, vertical outliers;
%%           3, leverage points;
%% Output: A, b

function [A,b] = Generate_A_b_Conf_large(m,n,flag_cm)
rng default;

a = 0.5;
outnum = 0.1*m; % the number of outliers in observations
%------------------------------ Generate matrix A  ------------------------
T = Generate_Toeplitz_matrix(a,n);
[R,idef] = chol(T);
if (idef)
    error('T is not positive definite');
end
A = randn(m,n)*R;
   
if (flag_cm == 3)
    A_tmp = 50 + randn(outnum,n); 
    outpos = randperm(m,outnum);
    A(outpos,:) = A_tmp;    
end

%-------------------------- Generate x_true -------------------------------
x_true = zeros(n,1);
x_true([1,3,7,10]) = 1.5;
x_true([2,9]) = 0.5;
x_true([4,6,8,11]) = 1;
Ax_ture = A*x_true;
%------------------------- Generate epsilon -------------------------------
epsilon = 0.5*randn(m,1);
if (flag_cm == 2) || (flag_cm == 3)
    epsilon_tmp = 20 + 0.5*randn(outnum,1);
    outpos = randperm(m,outnum);
    epsilon(outpos) = epsilon_tmp;
end
%------------------------------- Generate b -------------------------------
b = Ax_ture + epsilon;
%eval(['save A_b_',num2str(m),'_',num2str(n),'_cm',num2str(flag_cm),'.mat  A b']);

