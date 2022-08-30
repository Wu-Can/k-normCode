%%********************************************************************************************************
%% SIRPN_CD:
%% A smoothing method based on the inexact regularized proximal Newton
%% method whose subproblems are solved by the coordinate descent algorithm
%% for solving the convex CVaR-based sparse linear regression:
%% 
%%  minimize_{x in R^n, t in R} sum^n_{i=1}max(|b_i-A_i*x|-lambda*t,0) + lambda*k||(x,t)||_1
%%
%% where A is an n*d design matrix whose each row is A_i, b is the response vector in R^n,
%% lambda is the given positive parameter, k is the given positive integer between 1 and n.
%% [obj,x,t,runhist,info] = SIRPN_CD(Ainput,b,OPTIONS)
%% Input: 
%% Ainput, b
%% OPTIONS.epsilon = parameter in CHKS smoothing function;
%% OPTIONS.tol = the accuracy tolerance for solving the problem;
%% OPTIONS.tolsub = the accuracy tolerance for solving subproblems by IRPN
%% OPTIONS.obj_opt = the optimal value of the problem;
%% OPTIONS.maxiter = the maximum number of outer iteration;
%% OPTIONS.maxitersub = the maximum number of inner iteration;
%% OPTIONS.maxtime = the maximum time of the S-IRPN;
%% OPTIONS.mu = the parameter in Hessian for the subproblems of IRPN;
%% Output:
%% obj = objective value of the problem;
%% (x,t) = the variable of the problem;
%% runhist = a structure containing the history of the run;
%% info.iter = the total number of outer iteration;
%% info.itersub = the total number of inner iteration;
%% info.relobj = relative residual between obj and obj_opt;
%% info.time =  total time;
%% info.xnnz = the number of nonzeros in x;
%% SIRPN_CD:
%% Copyright (c) 2022 by
%% Can Wu, Ying Cui, Donghui Li, Defeng Sun
%%******************************************************************************************************

function [obj,x,t,runhist,info] = SIRPN_CD(Ainput,b,OPTIONS)

%%
%% Input parameters
%%
maxiter = 10;
maxitersub = 200;
maxtime = 7200;
tol = 1e-3;
tolsub = 1e-2;
tol_sparse = 0.14;
flag_sparse = 0;
UCI = 0;

if isfield(OPTIONS,'kk'), kk = OPTIONS.kk; end
if isfield(OPTIONS,'n'), n = OPTIONS.n; end
if isfield(OPTIONS,'d'), d = OPTIONS.d; end
if isfield(OPTIONS,'epsilon'), epsilon = OPTIONS.epsilon; end
if isfield(OPTIONS,'tol'), tol = OPTIONS.tol; end
if isfield(OPTIONS,'tolsub'), tolsub = OPTIONS.tolsub; end
if isfield(OPTIONS,'obj_opt'), obj_opt = OPTIONS.obj_opt; end
if isfield(OPTIONS,'maxiter'), maxiter = OPTIONS.maxiter; end
if isfield(OPTIONS,'maxitersub'), maxitersub = OPTIONS.maxitersub; end
if isfield(OPTIONS,'maxtime'), maxtime = OPTIONS.maxtime; end
if isfield(OPTIONS,'mu'), mu = OPTIONS.mu; end
if isfield(OPTIONS,'lambda'), lambda = OPTIONS.lambda;
else, lambda = lamc*max(abs(b'*Ainput)); end
if isfield(OPTIONS,'UCI'), UCI = OPTIONS.UCI; end

par.kk = kk;
par.n = n;
par.d = d;
par.lambda = lambda;
par.epsilon = epsilon;
par.mu = mu;

printyes = 1;
cntATmap = 0;
%%
%% Amap and ATmap
%%
if isstruct(Ainput)
    if isfield(Ainput,'A'); A = Ainput.A; end
    if isfield(Ainput,'ATmap'); ATmap = Ainput.ATmap; end
else
    A = Ainput;
    ATmap = @(y) (y'*A)'; %ATmap = @(x) mexAx(A,x,1);
end

sparseA = (n*d - nnz(A))/n*d;
if sparseA > tol_sparse
    A = sparse(A);
    flag_sparse = 1;
end
AT2 = (A.^2)';
%%
%% Generate initial point
%%
tstart = clock;
tstart_cpu = cputime;
x0 = zeros(d,1);
t0 = 0;
%%
%% print the initial information
%%
x_new = x0; t_new = t0; xt_new = [x_new;t_new];
b_Ax_new = b;
Ax_new = zeros(n,1);

obj_F = sum(max(abs(b_Ax_new)-par.lambda*t_new,0)) + (par.lambda*par.kk)*norm(xt_new,1);
relobj = abs(obj_F - obj_opt)/(1+abs(obj_opt));
ttime = etime(clock,tstart);

if printyes
    fprintf('\n *******************************************************');
    fprintf('****************************************************');
    fprintf('\n \t\t   SIRPN  for solving CVaR-based sparse linear regression ');
    fprintf('with k = %6.3f and lambda = %3.2e', par.kk,par.lambda);
    fprintf('\n ************************************************************');
    fprintf('***********************************************\n');
    fprintf('\n problem size: n = %3.0f, d = %3.0f', n, d);
    fprintf('\n obj_opt = %9.8e, obj_initial = %9.8f, relobj = %3.2e',obj_opt,obj_F,relobj);
    fprintf('\n epsilon0 = %3.2e', epsilon);
    fprintf('\n ---------------------------------------------------');
    fprintf('-------------------------');
    fprintf('\n  iters  [r_k       stoptol     iterCD]  [step   itersp]    time   [obj_Feps    del_Feps]    mu    |   resf');
    fprintf('      relobj');
end
%%
%% Main loop of SIRPN_CD
%%
for iter = 1:maxiter
    par.epsilon = 0.1*par.epsilon;
    %% Stopping criterion of SIRPN_CD
    if (relobj < tol) || (ttime > maxtime)
        if (relobj < tol)
            msg = 'SIRPN converged';
            info.termcode = 1;
        else
            msg = ' maximum time reached';
            info.termcode = 3;
        end
        break;
    end
     
    if iter < 3
       par.maxiterCD = 1000;
    elseif iter < 5
        par.maxiterCD = 2000;
    elseif iter < 6
        par.maxiterCD = 3000;
    else
        par.maxiterCD = 4000;
    end
    
    if iter == 1
        tolsub = 1e-1;
    else
        tolsub = 1e-2;
    end
    
    %%
    %% IRPN for inner subproblem
    %%
    x_ssub = x_new; t_ssub = t_new; xt_ssub = [x_ssub;t_ssub];
    parsub = par; break_inner = 0;
    b_Ax_ssub = b_Ax_new; Ax_ssub = Ax_new;
    numCD.iterCD = [];
    runhistsub.relobj = [];
    
    for itersub = 1:maxitersub
        phi_lamt_sub = phi_eps_fun(b_Ax_ssub,parsub.epsilon,0)-parsub.lambda*t_ssub;
        feps_ssub = sum(psi_eps_fun(phi_lamt_sub,parsub.epsilon,0));
        
        grad_psi = psi_eps_fun(phi_lamt_sub,parsub.epsilon,1);
        grad_phi = phi_eps_fun(b_Ax_ssub,parsub.epsilon,1);
        grad_f = -[ATmap(grad_psi.*grad_phi);parsub.lambda*sum(grad_psi)]; cntATmap = cntATmap + 1;
        resf = norm(grad_f + Proj_inf(xt_ssub - grad_f,parsub.lambda*parsub.kk));
  
        %% Stopping criterion of the inner subproblem 
        if ((resf < tolsub) || (ttime > maxtime)) && (itersub > 1)
            break_inner = 1;
            break;
        end
        
        %% Solving the subproblem minimize q_k(x,t) by CD
        hess_psi = psi_eps_fun(phi_lamt_sub,parsub.epsilon,2);
        hess_phi = phi_eps_fun(b_Ax_ssub,parsub.epsilon,2);
        beta_vec = hess_psi.*grad_phi;
        D_vec = beta_vec.*grad_phi + grad_psi.*hess_phi;
        D = sparse(1:parsub.n,1:parsub.n,D_vec); 
        lamATbeta = parsub.lambda*ATmap(beta_vec); cntATmap = cntATmap + 1;
        lam2hpsi = parsub.lambda^2*sum(hess_psi);
        DA = D*A;
        
        parsub.rho = 0.5; 
        parsub.mu = min(max(parsub.mu*max(resf,resf^(parsub.rho)),1e-12),1);
        Hxt_ssub = [((Ax_ssub)'*(DA))' + lamATbeta*t_ssub; lamATbeta'*x_ssub + lam2hpsi*t_ssub] + parsub.mu*xt_ssub;
        
        %% Coordinate descent algorithm for solving subproblems of IRPN
        gradf_Hess = grad_f - Hxt_ssub;
        ATDA_vec = AT2*D_vec;
        parsub.eta = 0.5;
        tol_CD = parsub.eta*min([resf,resf^(1 + parsub.rho)]);
        
        [xt_ssnew,Ax_ssnew,parsub] = ...
            CD(A,xt_ssub,tol_CD,resf,gradf_Hess,beta_vec,DA,lamATbeta,Ax_ssub,lam2hpsi,ATDA_vec,flag_sparse,parsub);
        numCD.iterCD(itersub) = parsub.iterCD;
        
        %% Backtracking line search
        dir = xt_ssnew - xt_ssub;
        g_ssub = (parsub.lambda*parsub.kk)*norm(xt_ssub,1);
        obj_Feps_ssub = feps_ssub + g_ssub;
        parsub.theta = 0.25; 
        parsub.beta = 0.25; 
        
        del_beta = 1; maxii = 20;
        for ii = 1:maxii
            del_xt_tmp = del_beta*dir;
            xt_tmp = xt_ssub + del_xt_tmp;
            x_tmp = xt_tmp(1:parsub.d);
            t_tmp = xt_tmp(end);
            
            % Compute the value of functions ell and F_eps at xt_tmp
            g_tmp = (parsub.lambda*parsub.kk)*norm(xt_tmp,1);
            ell_tmp = feps_ssub + grad_f'*del_xt_tmp + g_tmp;
            del_ell_tmp = parsub.theta*(obj_Feps_ssub - ell_tmp);
            
            Ax_tmp = (1 - del_beta)*Ax_ssub + del_beta*Ax_ssnew;
            b_Ax_tmp = b - Ax_tmp; 
            phi_lamt_tmp = phi_eps_fun(b_Ax_tmp,parsub.epsilon,0)-parsub.lambda*t_tmp;
            obj_Feps_tmp = sum(psi_eps_fun(phi_lamt_tmp,parsub.epsilon,0)) + g_tmp;
            del_Feps_tmp = obj_Feps_ssub - obj_Feps_tmp;
            
            if (del_ell_tmp <= del_Feps_tmp) || (ii == maxii)
                x_ssub = x_tmp;
                t_ssub = t_tmp;
                xt_ssub = xt_tmp;
                b_Ax_ssub = b_Ax_tmp;
                Ax_ssub = Ax_tmp;
                g_ssub = g_tmp;
                if ii == maxii
                    fprintf('* ');
                end
                break;
            end
             del_beta = del_beta*parsub.beta; % step length
        end
        
        obj_F = sum(max(abs(b_Ax_ssub)-parsub.lambda*t_ssub,0)) + g_ssub;
        relobjsub = abs(obj_F - obj_opt)/(1+abs(obj_opt));       
        if (relobjsub <= tol)
            msg = 'SIRPN converged';
            info.termcode = 1;
            break;
        end
        
        %%  Check stagnation for IRPN
        runhistsub.obj_F(itersub) = obj_F;
        if (itersub > 15) && all(runhistsub.obj_F(itersub-2:itersub) > runhistsub.obj_F(itersub-3:itersub-1))
            if (printyes);  fprintf('#'); end
            break;
        end
        
        if (itersub > 20)
            ratio = runhistsub.obj_F(itersub-9:itersub)./runhistsub.obj_F(itersub-10:itersub-1);
            if ((UCI == 1)&&(min(ratio) > 0.99999)) || ((UCI == 0)&&(min(ratio) > 0.99998))
                if (printyes);  fprintf('##'); end
                break;
            end
        end
        
        
        %% Print results of IRPN
        ttime = etime(clock,tstart);
        if (printyes)
            fprintf('\n [%4.0d]  [%3.2e  %3.2e     %5.0d]  [%3.2e   %2.0d]  %5.1f    [%6.4e  %3.2e] %3.2e | %3.2e',...
                itersub,parsub.resf_CD,tol_CD,parsub.iterCD,del_beta,ii,ttime,obj_Feps_ssub,del_Feps_tmp,parsub.mu,resf);
            fprintf('  %3.2e', relobjsub);
        end
        
    end
    x_new = x_ssub; t_new = t_ssub; xt_new = xt_ssub;
    b_Ax_new = b_Ax_ssub; Ax_new = Ax_ssub;
    obj_F = sum(max(abs(b_Ax_new)-par.lambda*t_new,0)) + (par.lambda*par.kk)*norm(xt_new,1);
    relobj = abs(obj_F - obj_opt)/(1+abs(obj_opt));
    
    ttime = etime(clock,tstart);
    if break_inner == 1
        runhist.itersub(iter) = itersub - 1;
    else
        runhist.itersub(iter) = itersub;
    end
    if (printyes)
        fprintf('\n iter = %5.0d, eps = %3.2e, relobj = %3.2e, obj_F = %6.6e',...
            iter,par.epsilon,relobj,obj_F);
    end
    
    runhist.obj_F(iter) = obj_F;
    runhist.stepiter(iter) = ii;
    runhist.relobj(iter) = relobj;
    runhist.iterCD(iter) = sum(numCD.iterCD);
end

x = x_new; t = t_new;
obj = obj_F;
if (iter == maxiter)
    msg = ' maximum iteration reached ';
    info.termcode = 2;
else
    iter = iter - 1;
end
%----------------------the number of nonzeros in x-------------------------
sortx = sort(abs(x),'descend');
normx1 = 0.999*norm(x,1);
tmpidex = find(cumsum(sortx) > normx1);
if isempty(tmpidex)
    nnzeros_x = 0;
else
    nnzeros_x = tmpidex(1);
end
%--------------------------------------------------------------------------
info.minx = min(x); info.maxx = max(x);
info.numiterCD = sum(runhist.iterCD);
info.cntATmap = cntATmap;
itersub_total = sum(runhist.itersub);
ttime = etime(clock,tstart);
ttime_cpu = cputime - tstart_cpu;
%%
%% Print Results
%%
if (printyes)
    if ~isempty(msg); fprintf('\n %s',msg); end
    fprintf('\n--------------------------------------------------------------');
    fprintf('------------------');
    fprintf('\n  number iter = %2.0d',iter);
    fprintf('\n  number iter of inner problems for IRPN = %2.0d',itersub_total);
    fprintf('\n  IRPN per iter = %3.0f', itersub_total/iter);
    fprintf('\n  number iterations of CD = %5.0d',info.numiterCD);
    fprintf('\n  time = %3.2f',ttime);
    fprintf('\n  time per iter = %5.4f',ttime/iter);
    fprintf('\n  cputime = %3.2f', ttime_cpu);
    fprintf('\n  objective value = %9.8e',obj);
    fprintf('\n  number of nonzeros in x = %3.0f', nnzeros_x);
    fprintf('\n  min(x) = %3.2e, max(x) = %3.2e', info.minx, info.maxx);
    fprintf('\n  relobj = %3.2e', relobj);
    fprintf('\n ATmap cnt = %4.0f', info.cntATmap);
    fprintf('\n--------------------------------------------------------------');
    fprintf('------------------\n');
end
info.xnnz = nnzeros_x;
info.relobj = relobj;
info.iter = iter;
info.itersub = itersub_total;
info.time = ttime;
info.obj_F = obj_F;


%--------------------------------------------------------------------------
