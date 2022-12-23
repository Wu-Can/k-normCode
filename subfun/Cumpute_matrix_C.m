%%=========================================================================
%% Cumpute_matrix_C:
%% compute the matrix C such that C*C^T = AJ*AJ^T + W
%%
%% C = Cumpute_matrix_C(flag_case,options)
%%
%% Input:
%% flag_case = the positive number depending on the form of the
%%             matrix W
%% options.AJ = the matrix AJ
%% options.index_alph = the index set alpha
%% options.index_beta = the index set beta
%% options.Palpha = submatrix of P by only keeping all the rows in alpha
%% options.Pbeta = submatrix of P by only keeping all the rows in beta
%% options.Pgamma = submatrix of P by only keeping all the rows in gamma
%% options.m = the sample size
%% Output:
%% C = the matrix C
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the Section 3.3 in the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%=========================================================================
function C = Cumpute_matrix_C(flag_case,options)
AJ = options.AJ;
index_alpha = options.index_alpha;
index_beta = options.index_beta;
Palpha = options.Palpha;
Pbeta = options.Pbeta;
Pgamma = options.Pgamma;
m = options.m;

switch(flag_case)
    case 1
        C = AJ;
    case 2
        C = [AJ, (Palpha)'];
    case {3.1, 3.2}
        C = [AJ, (1/sqrt(index_alpha))*(ones(1,index_alpha)*Palpha)', (Pbeta)'];
    case {4.1, 4.2}
        C = [AJ, (Palpha)', (1/sqrt(index_beta))*(ones(1,index_beta)*Pbeta)', (Pgamma)'];
    case 5
        C = [AJ,speye(m)];
end


