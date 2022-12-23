%%=========================================================================
%% Cumpute_matrix_W_MM:
%% calculate an element in the Clarke generalized Jacobian 
%% of Prox_{(1/sigma)*||\cdot||_{k_1}}(\cdot)
%%
%% W =  Cumpute_matrix_W_MM(flag_case,options)
%%
%% Input:
%% options.index_alpha = the index set alpha
%% options.index_beta = the index set beta
%% options.index_gamma = the index set gamma
%% options.Palpha = submatrix of P by only keeping all the rows in alpha
%% options.Pbeta = submatrix of P by only keeping all the rows in beta
%% options.Pgamma = submatrix of P by only keeping all the rows in gamma
%% options.m = the sample size
%% Output:
%% W = the matrix W
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%=========================================================================
function W =  Cumpute_matrix_W_MM(flag_case,options)
index_alpha = options.index_alpha;
index_beta = options.index_beta;
index_gamma = options.index_gamma;
Palpha = options.Palpha;
Pbeta = options.Pbeta;
Pgamma = options.Pgamma;
m = options.m;

switch(flag_case)
    case 1
        W = sparse(m,m);
    case {2.1, 2.2}
        W = (Palpha)'*Palpha;
    case {3.1, 3.2}
        if index_gamma == 0
            eTP = ones(1,index_alpha+index_beta)*[Palpha;Pbeta];
            W = 1/(index_alpha+index_beta)*(eTP'*eTP);
        else
            eTP = ones(1,index_alpha+index_beta)*[Palpha;Pbeta];
            W = 1/(index_alpha+index_beta)*(eTP'*eTP) + Pgamma'*Pgamma;
        end
    case {4.1, 4.2}
        if index_gamma == 0
            eTP = ones(1,index_beta)*Pbeta;
            W = Palpha'*Palpha + (1/index_beta)*(eTP'*eTP);
        else
            eTP = ones(1,index_beta)*Pbeta;
            W = Palpha'*Palpha + (1/index_beta)*(eTP'*eTP) + Pgamma'*Pgamma;
        end
    case 5
        W = speye(m);
end

end


