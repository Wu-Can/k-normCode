%%*************************************************************************
%% MM_NPPA:
%% A majorization-minimization algorithm associated the semismooth 
%% Newton method based on the proximal point algorithm for solving 
%% the nonconvex truncated CVaR-based sparse linear regression:
%%
%% (P) minimize_{x in R^n} {||A*x - b||_(k1) - ||A*x - b||_(k2) 
%%                          + lambda ||x||_1}
%% 
%% where A in R^{m * n}, b in R^m and lambda in R_+ are given,
%% and n >= k1 > k2 > 0.
%%
%% [x,obj,info,runhist] = MM_NPPA(A,b,kk1,kk2,lambda,OPTIONS)
%% 
%% Input:
%% A, b, kk1(k1), kk2(k2), lambda
%% OPTIONS.Warm_starting = 1, warm starting with the N-ALM for a 
%%                            good initial point;
%%                         0, without any warm starting and origin 
%%                            as the initial point;
%% OPTIONS.tol = accuracy tolerance for solving the problem;
%% OPTIONS.tol_alm = accuracy tolerance for the warm-starting N-ALM;
%% OPTIONS.maxiterMM = maximum number of outer iteration of the 
%%                     MM_NPPA;
%% OPTIONS.maxtime = maximum time for the MM_NPPA;
%% OPTIONS.rho0 = initial value of the parameter rho in the MM;
%% OPTIONS.sigma0 = initial value of the parameter sigma 
%%                  (c in manuscript) in the MM;
%% OPTIONS.eta0 = initial value of the parameter eta in the PPA;
%% OPTIONS.rhoscale = coefficient in (0,1) for updating rho;
%% OPTIONS.sigmascale = coefficient in (0,1) for updating sigma;
%% OPTIONS.etascale = coefficient in [1,+inf) for updating eta;
%% OPTIONS.UCI = 1, using the UCI data; 0, using the random data;
%% OPTIONS.rho_iter = frequency for updating rho and sigma;
%% Output:
%% x = the output solution x of (P);
%% obj = the output objective value of (P);
%% info.iter = total number of outer MM iteration;
%% info.iterPPA = total number of inner PPA iteration;
%% info.iterSSN = total number of inner SSN iteration;
%% info.obj_gap = relative residual for two adjacent objectives;
%% info.time = total running time;
%% info.nnzeros_x = the number of nonzero entries for x;
%% runhist = a structure containing the history of the run;
%% MM_NPPA:
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the Section 5 in the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%*************************************************************************

function [x,obj,info,runhist] = MM_NPPA(A,b,kk1,kk2,lambda,OPTIONS)

Warm_starting = 1;
tol = 1e-6;
tol_alm = 1e-4;
maxiterMM = 200;
maxtime = 7200; 
rho0 = 1;
sigma0 = 1;
eta0 = 1;
rhoscale = 0.5; 
sigmascale = 0.5;
etascale = 2; 

UCI = 0;
printMM = 1; % print in MM;
printPPA = 0; % print in PPA
printSSN = 0; % print in SSN
rho_iter = 4;
tiny = 1e-10;
Test_a = 0; % 1, print the information of the subgradient of ||Ax - b||_(kk2) 
%%
%% Input parameters
%%
if isfield(OPTIONS,'Warm_starting'), Warm_starting = OPTIONS.Warm_starting; end
if isfield(OPTIONS,'tol'), tol = OPTIONS.tol; end
if isfield(OPTIONS,'tol_alm'), tol_alm = OPTIONS.tol_alm; end
if isfield(OPTIONS,'maxiterMM'), maxiterMM = OPTIONS.maxiterMM; end
if isfield(OPTIONS,'maxtime'), maxtime = OPTIONS.maxtime; end
if isfield(OPTIONS,'rho0'), rho0 = OPTIONS.rho0; end
if isfield(OPTIONS,'sigma0'), sigma0 = OPTIONS.sigma0; end
if isfield(OPTIONS,'eta0'), eta0 = OPTIONS.eta0; end
if isfield(OPTIONS,'rhoscale'), rhoscale = OPTIONS.rhoscale; end
if isfield(OPTIONS,'sigmascale'), sigmascale = OPTIONS.sigmascale; end
if isfield(OPTIONS,'etascale'), etascale = OPTIONS.etascale; end
if isfield(OPTIONS,'rho_iter'), rho_iter = OPTIONS.rho_iter; end
if isfield(OPTIONS,'m'), m = OPTIONS.m; else, m = length(b); end
if isfield(OPTIONS,'n'), n = OPTIONS.n; else, n = size(A,2); end
if isfield(OPTIONS,'x0'), x0 = OPTIONS.x0; else, x0 = zeros(n,1); end
if isfield(OPTIONS,'u0'), u0 = OPTIONS.u0; else, u0 = zeros(m,1); end
if isfield(OPTIONS,'lamc'), lamc = OPTIONS.lamc; end
if isfield(OPTIONS,'UCI'), UCI = OPTIONS.UCI; end

par.m = m;
par.n = n;
par.kk1 = kk1;
par.lambda = lambda;
par.printPPA = printPPA;
par.printSSN = printSSN;
par.etascale = etascale;
par.rho = rho0;
par.sigma = sigma0;
cntATmap = 0;
cntAmap = 0;
breakPMM = 0;
iterPPA = 0;
iterSSN = 0;
obj_gap = inf;
%%
%% Amap and ATmap
%%
Amap = @(x) mexAx(A,x,0);
ATmap = @(y) (y'*A)';
%%
%% Generate initial point
%%
tstart = clock;
tstart_cpu = cputime;
%% Warm starting by N-ALM
if Warm_starting == 1
    fprintf('\n *******************************************************');
    fprintf('**********************************************************');
    fprintf('\n \t\t   Warm starting by N-ALM for the truncated CVaR-based model');
    fprintf('with kk1 = %4.0f, kk2 = %4.0f, lambda = %3.2e',kk1,kk2,lambda);
    fprintf('\n ******************************************************');
    fprintf('***********************************************************\n');
    par_alm.tol = tol_alm;
    par_alm.m = m;
    par_alm.n = n;
    par_alm.lambda = lambda;
    
    par_alm.sigmascale = 1.3;
    par_alm.tauscale = 1.2;
    if UCI == 1
        par_alm.sigma = 1e-4/lamc;
        par_alm.tau = 1e-4/lamc;
    else
        par_alm.sigma = 20/lambda;
        par_alm.tau = 20/lambda;
    end
    
    par_alm.x0 = x0;
    par_alm.u0 = u0;
    
    [~,x_alm,~,u_alm,~,info_alm] = NALM(A,b,par_alm); %  relkkt < tol_alm
    x = x_alm; u = u_alm;
    par.rho = 0.01;
    par.sigma = 0.01;
    info.warm_time = info_alm.totletime;
else
    x = x0; u = u0;
end
normx1 = norm(x,1);
if normx1 == 0
    Ax = sparse(m,1);
else
    Ax = Amap(x);
    cntAmap = cntAmap + 1;
end

Axb = Ax-b; ATu = ATmap(u);
Axbsort = sort(abs(Axb),'descend');
primobj_old = sum(Axbsort(1:kk1))-sum(Axbsort(1:kk2))+lambda*normx1;

%%
%% print the initial information
%%
if printMM
    fprintf('\n *******************************************************');
    fprintf('**********************************************************');
    fprintf('\n \t\t   MM+N-PPA for solving the truncated CVaR-based model with');
    fprintf('kk1 = %4.0f,kk2 = %4.0f, lambda = %3.2e', kk1,kk2,lambda);
    fprintf('\n ******************************************************');
    fprintf('***********************************************************\n');
    fprintf('\n rho0 = %3.4f, sigma0 = %3.4f, eta0 = %3.4f', par.rho, par.sigma, eta0);
    fprintf('\n rhoscale = %3.2f, sigmascale = %3.2f, etascale = %3.2f',...
        rhoscale, sigmascale, etascale);
    fprintf('\n -------------------------------------------------------');
    fprintf('----------------------------------------------------------\n');
    fprintf('\n problem size: m = %3.0f, n = %3.0f',m, n);
    fprintf('\n primal obj: %5.8e',primobj_old);
    fprintf('\n -----------------------------------------------------------');
    fprintf('-----------------------------------------------------------------');
    fprintf('\n iter  |  err_ppa    tol_ppa  |   relkkt  |    rho       sigma      eta    |');
    fprintf(' time  | nnz |    objnew     objold    |  obj_gap  ');
    if Test_a == 1
        fprintf(' |  ell0  k  ell1 |    zk    | norm(a-aold) ');
    end
    
end

%%
%% ------------------------ Begin MM ------------------------------
%%
for iter = 1:maxiterMM
    
    %% choose one subgradient a from the subdifferential of norm(Ax-b,kk2)
    z = Axb;
    p = ones(m,1);
    p(z < -tiny) = -1; % avoid p(i) = 0 if z(i) = 0;
    [z_tmp,id] = sort(abs(z),'descend');
    P = sparse(1:m,id,p(id)); 
    
    z_tmpkk2 = z_tmp(kk2);
    index1 = ((z_tmp - z_tmpkk2) > tiny);
    k0 = sum(index1);
    index3 = ((z_tmp - z_tmpkk2) < -tiny);
    index2 = (abs(z_tmp - z_tmpkk2) <= tiny);
    k1 = k0 + sum(index2);
    
    a_tmp = ones(m,1);
    if z_tmpkk2 > tiny
        a_tmp(index3) = 0;
        if k1 > kk2
            a_tmp(kk2+1:k1) = 0;
        end
    else
        if kk2 < m
            a_tmp(kk2+1:m) = 0;
        end
    end
    a = P'*a_tmp;
    par.a = a;
   
    par.x0 = x;
    par.u0 = u;
    par.ATu = ATu;
    par.Ax = Ax;
    
    par.eta0 =  min(max(eta0,1/(115*par.rho)),100); 
    par.obj_gap_mm = obj_gap;
    if (iter > 1) && (Test_a == 1)
        difa_aold = norm(a-a_old);
    end
    a_old = a;
  
    %% compute [u, x, z] by the dual PPA
    [x,u,~,info_ppa,runhist_ppa] = dPPA_SSN(A,b,par);
    
    Axb = runhist_ppa.Ax_b;
    Ax = runhist_ppa.Ax;
    Axbsort = sort(abs(Axb),'descend');
    primobj_new = sum(Axbsort(1:kk1))-sum(Axbsort(1:kk2))+lambda*norm(x,1);
    def_primobj = primobj_new-primobj_old;
    obj_gap = abs(def_primobj)/(1+abs(primobj_old));
   
    %% stopping criterion
    if obj_gap  < tol
        msg = 'converged';
        breakPMM = 1;
    end
    ttime = etime(clock,tstart);
    
    %% Print the result of dual PPA+SSN
    if printMM == 1
        fprintf('\n  %3d  |  %3.2e   %3.2e |  %3.2e | %3.2e   %3.2e   %3.2e |',...
            iter, info_ppa.err_ppa, info_ppa.tol_ppa,...
            info_ppa.relkkt, par.rho, par.sigma, info_ppa.eta);
        fprintf(' %5.2f | %3.0f | %- 5.4e %- 5.4e |  %3.2e ',...
            ttime, info_ppa.nnzeros_x, primobj_new, primobj_old, obj_gap);

        if (Test_a == 1)
            fprintf('|  %3.0d  %3.0d  %3.0d | %3.2e |',k0,kk2,k1,z_tmpkk2);
            if iter > 1
                fprintf('%3.2e',difa_aold);
            end
        end
    end
    runhist.rho(iter) = par.rho;
    runhist.sigma(iter) = par.sigma;
    runhist.eta(iter) = info_ppa.eta;
    runhist.cultime(iter) = ttime;
    
    %% updata parameter   
    if mod(iter,rho_iter) == 0
        par.rho = max(rhoscale*par.rho,1e-6);
        par.sigma = max(sigmascale*par.sigma,1e-6);
    end
    
    %% record running history
    primobj_old = primobj_new;
    norma = sum(a_tmp);
    
    runhist.objPMM(iter) = primobj_new;
    runhist.norma(iter) = norma;
    runhist.a{iter} = a;
    runhist.z_tmpkk2_vec(iter) = z_tmpkk2;
    runhist.obj_gap(iter) = obj_gap;
    runhist.iterSSN(iter) = info_ppa.iterSSN;
    runhist.iterPPA(iter) = info_ppa.iterPPA;
    runhist.time(iter) = info_ppa.time;
    
    cntATmap = cntATmap + info_ppa.cntATmap;
    cntAmap = cntAmap + info_ppa.cntAmap;
    ATu = runhist_ppa.ATu;
    
    iterPPA = iterPPA + info_ppa.iterPPA;
    iterSSN = iterSSN + info_ppa.iterSSN;
    %%
    if  (iter == maxiterMM)
        msg = ' maximum iteration reached';
        breakPMM = 100;
    end
    if  (ttime > maxtime)
        msg = ' maximum time reached';
        breakPMM = 1000;
    end
    if def_primobj > 1e-10
        msg = ' objective is increasing';
        breakPMM = 10000;
    end
    
    if (breakPMM > 0) || (ttime > maxtime) || (iter == maxiterMM)
        fprintf('\n  breakyes = %3.1f, %s, obj_gap = %3.2e',breakPMM,msg,obj_gap);
        break;
    end
end

%%
%% ------------------------ End MM --------------------------------
%%
obj = primobj_new;

ttime = etime(clock,tstart);
ttime_cpu = cputime - tstart_cpu;
%----------------------the number of nonzeros in x-------------------------
[sortx,index_x] = sort(abs(x),'descend');
normx1 = 0.999*norm(x,1);
tmpidex = find(cumsum(sortx) > normx1);
if isempty(tmpidex)
    nnzeros_x = 0;
else
    nnzeros_x = tmpidex(1);
end
info.minx = min(x); info.maxx = max(x);
index_nnz = index_x(1:nnzeros_x);
if nnzeros_x == n
    index_0 = [];
else
    index_0 = index_x(nnzeros_x+1:n);
end
info.totle_cntATmap = cntATmap + 1;
info.totle_cntAmap = cntAmap;
%--------------------------------------------------------------------------
%%
%% Print Results
%%
if (printMM)
    if ~isempty(msg); fprintf('\n  %s',msg); end
    fprintf('\n ------------------------------------------------------------');
    fprintf('\n  number iterPMM = %2.0d',iter);
    fprintf('\n  number iterPPA = %2.0d',iterPPA);
    fprintf('\n  number iterSSN = %2.0d',iterSSN);
    if Warm_starting == 1
        fprintf('\n  warm starting time = %3.2f',info.warm_time);
    end
    fprintf('\n  time = %3.2f',ttime);
    fprintf('\n  time per iter = %5.4f',ttime/iter);
    fprintf('\n  cputime = %3.2f', ttime_cpu);
    fprintf('\n  primobj = %9.8e',obj);
    fprintf('\n  obj_gap = %3.2e',obj_gap);
    fprintf('\n  number of nonzeros in x = %3.0f', nnzeros_x);
    fprintf('\n  min(x) = %3.2e, max(x) = %3.2e', info.minx, info.maxx);
    fprintf('\n  Amap cnt = %4.0f, ATmap cnt = %4.0f', info.totle_cntAmap,info.totle_cntATmap);
    fprintf('\n ------------------------------------------------------------\n');
end

%%
%% record history
%%
info.obj = obj;
info.iter = iter;
info.time = ttime;
info.x = x;
%info.z = z0;
%info.u = u0;
info.obj_gap = obj_gap;
info.nnzeros_x = nnzeros_x;
info.index_nnz = index_nnz;
info.index_0 = index_0;
info.z_tmpkk2 = min(runhist.z_tmpkk2_vec);
info.iterPPA = iterPPA;
info.iterSSN = iterSSN;
end


