%%===================================================================================
%% Generatedash_Jacobi:
%% obtain the information of the generalized Jacobian matrix
%% of projection on 
%%          B_*:={ z | \|z\|_inf <= r, \|z\|_1 <= kk*r }
%%
%% [flag_case,r_tmp,options] = Generatedash_Jacobi(A,w1,w2,options)
%%
%% Input:
%% A =  the m*n-dimensional design matrix 
%% w1 = the vector omega_1 defined in Section 3.3 of the paper
%% w2 = the vector omega_2 defined in Section 3.3 of the paper
%% options.kk = the parameter kk on B_*
%% options.lambda = the parameter lambda in the model (9) of the paper
%% options.sigma = the parameter sigma in Section 3.3 of the paper
%% options.m = the sample size
%% options.n = the feature size
%% options.indexJ =the index set indexJ defined in Section 3.3 of the paper
%% options.iter = the number of iteration of Algorithm 1 of the paper
%% Output:
%% flag_case = the a positive number depending on the form of the matrix D
%% r_tmp = the sum of the cardinality of several index sets
%% options.r_indexJ = the cardinality of indexJ
%% options.AJ = the matrix AJ
%% options.indexJ = the new index set indexJ
%% options.same_indexJ = 1, indexJ has not changed since the last SSN iteration
%%                     = 0, otherwise
%% options.index_alpha = the index set alpha
%% options.index_beta = the index set beta
%% options.index_gamma = the index set gamma
%% options.Palpha = the submatrix of P by only keeping all the rows in alpha
%% options.Pbeta = the submatrix of P by only keeping all the rows in beta
%% options.Pgamma = the submatrix of P by only keeping all the rows in gamma
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the Section 3.3 in the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%==================================================================================
function [flag_case,r_tmp,options] = Generatedash_Jacobi(A,w1,w2,options)
kk = options.kk;
lambda = options.lambda;
sigma = options.sigma;
n = options.n;
m = options.m;
indexJ_old = options.indexJ;
iter = options.iter;
tol = 1e-12;
r1 = sigma*lambda;
r = sigma;
indexJ_new = (abs(w1) > r1 + tol);
r_indexJ = sum(indexJ_new);
index_alpha = -1;
index_beta = -1;
index_gamma = -1;

if r_indexJ == 0
    AJ = sparse(m,n);
elseif r_indexJ < n
    if (iter > 1) && (norm(indexJ_new-indexJ_old,1) < 1)
        AJ = options.AJ;
    else
        AJ = A(:,indexJ_new);
    end
else
    AJ = A;
end

p = sign(w2);
[x,id] = sort(abs(w2),'descend');
P = sparse(1:m,id,p(id)); % save time

c = kk*r;
b = r*ones(m,1);
xp = Proj_dknorm(x,c,b);
tol = 1e-12;

if norm(indexJ_new-indexJ_old,1) < 1
    same_indexJ = 1;
else
    same_indexJ = 0;
end

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
if (flag_case == 1) || (relnormx_xp <= tol)
    flag_case = 1;
    r_tmp = r_indexJ;
    options.index_alpha = index_alpha;options.index_beta = index_beta;options.index_gamma = index_gamma;
    options.Palpha = []; options.Pbeta = []; options.Pgamma = [];
    return;
else
    switch (flag_case)
        case 2
            alp2 = (abs(r-xp) < tol); beta2 = ((abs(r-xp) >= tol));
            const21 = sum(alp2); const22 = sum(beta2);
            Palp2 = P(alp2,:); Pbeta2 = P(beta2,:);
            r_tmp = r_indexJ+const21;
            options.index_alpha = const21;options.index_beta = const22;options.index_gamma = index_gamma;
            options.Palpha = Palp2; options.Pbeta = Pbeta2; options.Pgamma = [];
            return;
        case 3
            alp3 = (abs(xp) > tol); const31 = sum(alp3);
            bet3 = (abs(xp) <= tol); const32 = sum(bet3);
            if const32 == 0
                Palp3 = P;
                Pbet3 = [];
                flag_case = 3.1;
            else
                Palp3 = P(alp3,:);
                Pbet3 = P(bet3,:);
                flag_case = 3.2;
            end
            r_tmp = r_indexJ+1+const32;
            options.index_alpha = const31; options.index_beta = const32;options.index_gamma = index_gamma;
            options.Palpha = Palp3; options.Pbeta = Pbet3; options.Pgamma = [];
            return;
        otherwise
            alp4 = (abs(r-xp) < tol); const41 = sum(alp4);
            bet4 = (xp > tol & abs(r-xp) >= tol); const42 = sum(bet4);
            gam4 = (xp <= tol); const43 = sum(gam4);
            if (const42 > 0)
                Palp4 = P(alp4,:);
                Pbet4 = P(bet4,:);
                if (const43 == 0)
                    Pgam4 = [];
                    flag_case = 4.1;
                else
                    Pgam4 = P(gam4,:);
                    flag_case = 4.2;
                end
                r_tmp = r_indexJ + const41 + 1 + const43;
                options.index_alpha = const41; options.index_beta = const42; options.index_gamma = const43;
                options.Palpha = Palp4; options.Pbeta = Pbet4; options.Pgamma = Pgam4;
                return;
            else
                flag_case = 5;
                r_tmp = r_indexJ + m;
                options.index_alpha = index_alpha;options.index_beta = index_beta;options.index_gamma = index_gamma;
                options.Palpha = []; options.Pbeta = []; options.Pgamma = [];
                return;
            end
    end
end




