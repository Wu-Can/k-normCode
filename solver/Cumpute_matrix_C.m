function C = Cumpute_matrix_C(flag_case,options)
%% This function is to compute the matrix C such that C*C^T = AJ*AJ^T + W
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


