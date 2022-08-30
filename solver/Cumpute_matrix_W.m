function W = Cumpute_matrix_W(flag_case,options)
%% This function is to compute the matrix W, which is an element of 
%% the Clarke generalized Jacobian of Prox_{sigma*||\cdot||_{k}}(\cdot)
%AJ = options.AJ;
index_alpha = options.index_alpha;
index_beta = options.index_beta;
Palpha = options.Palpha;
Pbeta = options.Pbeta;
Pgamma = options.Pgamma;
m = options.m;

switch(flag_case)
    case 1
        W = sparse(m,m);
    case 2
        W = (Palpha)'*Palpha;
    case 3.1
        eTPalp = ones(1,index_alpha)*Palpha;
        W = (1/index_alpha)*(eTPalp'*eTPalp);
    case 3.2
        eTPalp = ones(1,index_alpha)*Palpha;
        W = (1/index_alpha)*(eTPalp'*eTPalp) + (Pbeta)'*Pbeta;
    case 4.1
        eTPbeta = ones(1,index_beta)*Pbeta;
        W = Palpha'*Palpha + (1/index_beta)*(eTPbeta'*eTPbeta);
    case 4.2
        eTPbeta = ones(1,index_beta)*Pbeta;
        W = Palpha'*Palpha + (1/index_beta)*(eTPbeta'*eTPbeta) + Pgamma'*Pgamma;
    case 5
        W = speye(m);
end


