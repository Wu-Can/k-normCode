%%=========================================================================
%% matvec_N_ALM:
%% compute matrix*vector: Hess*d
%% where Hess = sigma*C*CT + (tau/sigma)*I_n
%%
%% Hessianv = matvec_N_ALM(d,C,options)
%%
%% Input:
%% d = a vecor
%% C = a part of matrix A (C*C^T = AJ*AJ^T + W)
%% options.sigma = the parameter sigma in Algorithm 2 of the paper
%% options.tau = the parameter tau in Algorithm 2 of the paper
%% Output:
%% Hessianv = the value of Hess*d
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the paper:
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%=========================================================================
function Hessianv = matvec_N_ALM(d,C,options)
sigma = options.sigma;
tau = options.tau;
Cmap = @(x) mexAx(C,x,0);
CTmap = @(x) mexAx(C,x,1);
Hessianv = sigma*Cmap(CTmap(d))+(tau/sigma)*d;
end