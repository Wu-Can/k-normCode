%%=========================================================================
%% CD:
%% coordinate descent algorithm (CD) for solving subproblems 
%% of S-IRPN
%%
%% [xt,Ax,options] = CD(A,xt_old,tol,resf,gradf_Hess,beta_vec,DA,...
%%           lamATbeta,Ax_new,lam2hpsi,ATDA_vec,flag_sparse,options)
%%
%% Input:
%% A = design matrix 
%% xt_old = current point (x, t) in the S-IRPN
%% tol = tolerance of stopping criterion of the CD
%% resf = the KKT residual of the subproblem of the S-IRPN at xt_old
%% gradf_Hess = difference between the gradient of f
%%              and Hessain*xt_old
%% beta_vec = hess_psi.*grad_phi
%% DA = D*A
%% lamATbeta = lambda*(AT*beta_vec)
%% Ax_new = A*x_old
%% lam2hpsi = lambda^2*sum(hess_psi)
%% ATDA_vec = AT2*D_vec with AT2 = (A.^2)' and D_vec = diag(D)
%% flag_sparse = 1, if A is sparse
%%             = 0, otherwise
%% options.maxiterCD = the maximum number of iterations of the CD
%% options.lambda = parameter lambda in (A8) of the supplementary 
%%                  materials
%% options.kk = parameter k in (A8) of the supplementary materials
%% options.d = feature size
%% options.mu = parameter mu in the S-IRPN
%% Output:
%% xt = the latest point (x, t) in the S-IRPN
%% Ax = A*x
%% options.iterCD = the number of iterations of the CD
%% options.resf_CD = the KKT residual of the subproblem of the
%%                   S-IRPN at xt
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the Appendix E in the supplementary 
%% materials of the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%=========================================================================
function [xt,Ax,options] = CD(A,xt_old,tol,resf,gradf_Hess,beta_vec,DA,...
    lamATbeta, Ax_new,lam2hpsi,ATDA_vec,flag_sparse,options)
maxiter = 1000;
if isfield(options,'maxiterCD'), maxiter = options.maxiterCD; end
lambda = options.lambda;
kk = options.kk;
lamkk = lambda*kk;
d = options.d;
mu = options.mu;
resf_sub = resf;
xt_new = xt_old;
Hess_vec = [ATDA_vec;lam2hpsi] + mu;

for iter = 1:maxiter
    if (mod(iter-1,10) == 0) && (resf_sub < tol)
        xt = xt_new;
        Ax = Ax_new;
        options.iterCD = iter - 1;
        options.resf_CD = resf_sub;
        break;
    end
    %% Compute xt_new
    if flag_sparse == 1
        [xt_new,Ax_new] = OnceCD(A,xt_new,ATDA_vec,DA,Ax_new,lamATbeta,gradf_Hess,lamkk,Hess_vec);
    else
        for jj = 1:d
            alpha = xt_new(jj)*ATDA_vec(jj) - dot(DA(:,jj),Ax_new) - xt_new(end)*lamATbeta(jj) - gradf_Hess(jj);
            xt_tmp = (alpha - max(-lamkk,min(lamkk,alpha)))/Hess_vec(jj);
            del_xt = xt_tmp - xt_new(jj);
            if abs(del_xt) > eps
                Ax_new = Ax_new + del_xt*A(:,jj);
            end
            xt_new(jj) = xt_tmp;
        end
    end
    lambetaTAx = lambda*(beta_vec'*Ax_new);
    alpha = -lambetaTAx - gradf_Hess(end);
    xt_tmp = (alpha - max(-lamkk,min(lamkk,alpha)))/Hess_vec(end);
    xt_new(end) = xt_tmp;
    
    if (mod(iter,10) == 0)
        Hxt = [((Ax_new)'*DA)'+ xt_tmp*lamATbeta; lambetaTAx + xt_tmp*lam2hpsi] + mu*xt_new + gradf_Hess;
        resf_sub = norm(Hxt + Proj_inf(xt_new - Hxt,lamkk));
    end
end
if iter == maxiter
    xt = xt_new;
    Ax = Ax_new;
    options.iterCD = iter;
    options.resf_CD = resf_sub;
end