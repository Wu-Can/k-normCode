%%=========================================================================
%% Generatedash_Jacobi_MM:
%% obtain some information of the generalized Jacobian matrix of 
%% projection on
%%        B_*={ z | \|z\|_inf<=r, \|z\|_1 <= kk*r }
%%
%% [flag_case,r_tmp,options] = Generatedash_Jacobi_MM(A,w1,w2,options)
%%
%% Input:
%% A = a m*n-dimensional design matrix 
%% w1 = a n-dimensional vector 
%% w2 = a m-dimensional vector 
%% options.lambda = the parameter lambda in the model (23) of the paper
%% options.rho = the parameter rho in (24) of the paper
%% options.sigma = the parameter c in (24) of the paper
%% options.m = the sample size
%% options.n = the feature size
%% options.indexJ = the index set indexJ
%% options.iter = the number of iteration of the PPA
%% options.AJ = the matrix AJ
%% Output
%% options.r_indexJ = the cardinality of indexJ
%% options.AJ = the new matrix AJ
%% options.indexJ = index set indexJ
%% options.same_indexJ = 1, indexJ has not changed since the last 
%%                          SSN iteration
%%                     = 0, otherwise
%% options.index_alpha = the index set alpha
%% options.index_beta = the index set beta
%% options.index_gamma = the index set gamma
%% options.gamISzero = 1, if const_gam = 0
%%                   = 0, otherwise
%% options.Palpha = submatrix of P by only keeping all the rows in alpha
%% options.Pbeta = submatrix of P by only keeping all the rows in beta
%% options.Pgamma = submatrix of P by only keeping all the rows in gamma
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%=========================================================================
function [flag_case,r_tmp,options] = Generatedash_Jacobi_MM(A,w1,w2,options)
kk = options.kk;
lambda = options.lambda;
rho = options.rho;
sigma = options.sigma;

n = options.n;
m = options.m;
indexJ_old = options.indexJ;
iter = options.iter;
tol = 1e-12; %tol = eps;
same_indexJ = 0;
r1 = lambda/rho;
r = 1/sigma;
theta2 = 1;
gamISzero = 0;

indexJ_new = (abs(w1) > r1 + tol);

r_indexJ = sum(indexJ_new);
r_tmp = -1;
if r_indexJ == 0
    AJ = sparse(m,n);
elseif r_indexJ < n
    if (iter > 1) && (norm(indexJ_new-indexJ_old,1) < 1)
        AJ = options.AJ;
        same_indexJ = 1;
    else
        AJ = A(:,indexJ_new);
    end
else
    AJ = A;
end

p = sign(w2);
[x,id] = sort(abs(w2),'descend');
P = sparse(1:m,id,p(id));

c = kk*r;
b = r*ones(m,1);
xp = Proj_dknorm(x,c,b);

options.r_indexJ = r_indexJ;
options.AJ = AJ;
options.indexJ = indexJ_new;
options.same_indexJ = same_indexJ;
normxp_inf = xp(1);
normxp_1 = sum(xp);

relr_xp_inf = (r-normxp_inf)/max(1,r);
relc_xp_1 = (c-normxp_1)/max(1,c);
relnormx_xp = sum(x-xp)/max(1,normxp_1);

if relr_xp_inf > tol
    if relc_xp_1 > tol
        flag_case = 1;
    else
        flag_case = 3;
    end
else
    if relc_xp_1 > tol
        flag_case = 2;
    else
        flag_case = 4;
    end
end
%--------------------------------------------------------------------------
dex_alpha = (abs(r-xp) < tol); const_alp= sum(dex_alpha);
dex_beta = (xp > tol & abs(r-xp) >= tol); const_bet = sum(dex_beta);
dex_gamma = (xp <= tol); const_gam = sum(dex_gamma);

Palpha = P(dex_alpha,:); Pbeta = P(dex_beta,:); Pgamma = P(dex_gamma,:);
if const_gam == 0
    gamISzero = 1;
end

options.index_alpha = const_alp;
options.index_beta = const_bet;
options.index_gamma = const_gam;
options.gamISzero = gamISzero;
options.Palpha = Palpha; options.Pbeta = Pbeta; options.Pgamma = Pgamma;
%--------------------------------------------------------------------------

if (relnormx_xp <= tol)
    flag_case = 1;
    r_tmp = r_indexJ;
    return;
else
    switch(flag_case)
        case 2
            r_tmp = r_indexJ + const_alp;
            if const_alp >= theta2*(const_bet+const_gam)
                flag_case = 2.1;
            else
                flag_case = 2.2;
            end
            return;
        case 3
            if (const_alp+const_bet)<=theta2*const_gam
                flag_case = 3.1;
            else
                flag_case = 3.2;
            end
            r_tmp = r_indexJ + 1 + const_gam;
            return;
        case 4
            if (const_bet > 0)
                if const_bet < theta2*(const_alp + const_gam)
                    flag_case = 4.1;
                else
                    flag_case = 4.2;
                end
                r_tmp = r_indexJ + const_alp + 1 + const_gam;
                return;
            else
                flag_case = 5;
                r_tmp = r_indexJ + m;
                return;
            end
    end
end





