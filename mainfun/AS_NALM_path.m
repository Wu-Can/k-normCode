%%*************************************************************************
%% AS_NALM_path: an adaptive sieving strategy associated N-ALM for 
%% computing a solution path of the convex CVaR-based sparse linear 
%% regression:
%%
%% (P)  minimize_{x in R^n} {||Ax - b||_(k) + lambda ||x||_1} 
%% where A in R^{m * n}, b in R^m and lambda > 0 are given.
%%
%% [Obj,X,U,runhistAS,infoAS] = AS_NALM_path(Ainput,b,OPTIONS)
%% Input:
%% Ainput,b;
%% OPTIONS.tol = accuracy tolerance for solving the problem;
%% OPTIONS.maxitersub = maximum number of inner iteration;
%% OPTIONS.lambda = a sequence of parameters lambda in (P) with 
%%                  length lenlam;
%% OPTIONS.kk = parameter k in (P);
%% OPTIONS.sigma = a sequence of the initial values for sigma in N-ALM;
%% OPTIONS.tau = a sequence of the initial values for tau in N-ALM;
%% OPTIONS.sigmascale = the positive number for updating sigma in N-ALM;
%% OPTIONS.tauscale = the positive number for updating tau in N-ALM;
%% Output:
%% Obj = a lenlam*2 matrix whose i-th row is the output primal and 
%%       dual objectives associated to lambda(i);
%% X = an n*lenlam matrix whose j-th column is the output primal 
%%     solution x associated to lambda(j);
%% U = a m*lenlam matrix whose j-th column is the output dual 
%%     solution u associated to lambda(j);
%% runhistAS.relkkt = a sequence of relative KKT residuals;
%% runhistAS.relgap = a sequence of relative duality gaps;
%% runhistAS.primfeas = a sequence of primal infeasibilities;
%% runhistAS.dualfeas = a sequence of dual infeasibilities;
%% runhistAS.iterAS = a sequence of the number of inner iterations 
%%                    for each outer AS iteration;
%% runhistAS.iterPAL = a sequence of the number of ALM iterations in 
%%                     N-ALM for each outer AS iteration;
%% runhistAS.iterSSN = a sequence of the number of SSN iterations in 
%%                     N-ALM for each outer AS iteration;
%% runhistAS.time = a sequence of the running time for each outer AS
%%                  iteration;
%% runhistAS.ttime_path = a sequence of the cumulative time;
%% runhistAS.xnnz_path = a sequence of the true cardinality of the 
%%                       solution
%% runhistAS.n_mean = a sequence of the average number of selected 
%%                    active features;
%% infoAS = a structure containing the history of the last AS 
%%          iteration; 
%% N-ALM:
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the Section 4 in the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%*************************************************************************
function [Obj,X,U,runhistAS,infoAS] = AS_NALM_path(Ainput,b,OPTIONS)
%%
%% Parameter setting
%%
maxitersub = 100;
tol = 1e-6;
lamc = 1e-3;

printyes = 1; % print in AS
print_yes = 0; % print in the outer ALM for N-ALM
print_level = 0; % print in the inner SSN for N-ALM
flag_tol = 2; % adopt relkkt < tol as the stopping criterion of the N-ALM
Test_stepALM = 1; % print the muber of iterations of N-ALM

if isfield(OPTIONS,'tol'), tol = OPTIONS.tol; end
if isfield(OPTIONS,'kk'), kk = OPTIONS.kk; end
if isfield(OPTIONS,'m'), m = OPTIONS.m; end
if isfield(OPTIONS,'n'), n = OPTIONS.n; end
if isfield(OPTIONS,'lambda'), lambda = OPTIONS.lambda;else, lambda = lamc*max(abs(b'*Ainput)); end
if isfield(OPTIONS,'sigma'), sigma = OPTIONS.sigma; end
if isfield(OPTIONS,'tau'), tau = OPTIONS.tau; end
if isfield(OPTIONS,'maxitersub'), maxitersub = OPTIONS.maxitersub; end
if isfield(OPTIONS,'flag_tol'), flag_tol = OPTIONS.flag_tol; end
if isfield(OPTIONS,'sigmascale'), sigmascale = OPTIONS.sigmascale; end
if isfield(OPTIONS,'tauscale'), tauscale = OPTIONS.tauscale; end

par.kk = kk;
par.m = m;
par.n = n;
par.tol = tol;
par.sigmascale = sigmascale;
par.tauscale = tauscale;
par.printyes = print_yes; 
par.printlevel = print_level; 
par.flag_tol = flag_tol;

msg = [];
%%
%% Amap and ATmap
%%
Amap = @(x) mexAx(Ainput,x,0);
ATmap = @(y) (y'*Ainput)';
%%
%% Print initial point of AS
%%
if printyes
    fprintf('\n *******************************************************');
    fprintf('******************************************');
    fprintf('\n \t\t   Adaptive sieving for solving k_norm with k=%6.3f ', par.kk);
    fprintf('\n ******************************************************');
    fprintf('*******************************************\n');
    fprintf('\n problem size: m = %3.0f, n = %3.0f',m, n);
    fprintf('\n ---------------------------------------------------');
    %fprintf('\n  iter|   Res    |priminfeas  dualinfeas   relgap |    pobj       dobj      |');
    %fprintf(' time |   sigma     tau   | res_kkt  | iterPAL   n      nnz');
end

%%
%% generate initial point from dADMM x0, z0, u0
%%
tstart = clock;
normb = norm(b);

% generate initial index set
s = abs(b'*Ainput)'./(sqrt(sum(Ainput.*Ainput)')*normb);
n_index = ceil(sqrt(n)/5);
s_new = sort(s,'descend');
s_num = s_new(n_index);
indexIb = (s >= s_num);
indexI = ~indexIb;

n_new = sum(indexIb);
x0 = zeros(n_new,1);
z0 = zeros(m,1);
u0 = z0;

par.x0 = x0; par.z0 = z0; par.u0 = u0;
lengthI = sum(indexI);
Apart = Ainput(:,indexIb);
maxiter = length(lambda);
par.n = n_new; normu = norm(par.u0);

Obj = zeros(maxiter,2);
X = zeros(n,maxiter);
U = zeros(m,maxiter);
Time = zeros(maxiter,1);
%%
%% Adaptive sieving strategy
%%
for iter = 1:maxiter
    par.lambda = lambda(iter);
    par.sigma = sigma(iter);
    par.tau = tau(iter);
    if normu == 0
        Res = 0;
    else
        ATu = ATmap(par.u0); normATu = norm(ATu);
        Res = norm(ATu + Proj_inf(x-ATu,par.lambda))/(1+norm(x)+normATu);
    end
    %par = OPTIONS;
    if (printyes)
        fprintf('\n Prob = %3.0d, lambda = %6.3f, Res = %3.2e, n = %6.0d', iter, par.lambda, Res, par.n);
        fprintf('\n ---------------------------------------------------');
        fprintf('\n  iter|   Res    |priminfeas  dualinfeas   relgap |    pobj       dobj      |');
        fprintf(' time |   sigma     tau   | res_kkt  ');  
        if Test_stepALM == 1
            fprintf('| numPAL  numSSN   n      nnz');
        end
    end
    par.tol = tol;

    [obj,x_part,z_old,u_new,~,info] = NALM(Apart,b,par);
      
    z_tmp = u_new + z_old;
    z_new = z_tmp - Proj_dknorm(z_tmp,par.kk,ones(m,1));
    ATu_part = info.ATu;
    x_tmp = x_part-ATu_part;
    x_partnew = x_tmp-Proj_inf(x_tmp,par.lambda);
    x_new = zeros(n,1);
    x_new(indexIb) = x_partnew;
    ATu = ATmap(u_new);
    Res = norm(ATu + Proj_inf(x_new-ATu,par.lambda))/(1+norm(x_new)+norm(ATu));
    
    if (printyes)
        fprintf('\n %5.1d| %3.2e |%3.2e    %3.2e   %- 3.2e| %- 5.4e %- 5.4e |',...
            0,Res,info.priminfeas_final,info.dualinfeas_final,info.res_gap_final,...
            obj(1),obj(2));
        fprintf('%5.1f | %3.2e %3.2e | %3.2e',info.totletime, info.sigma, info.tau,info.res_kkt_final);
        if Test_stepALM == 1
            fprintf(' | %5.1d  %5.1d  %5.1d  %5.1d',info.iter, info.numSSN, par.n, info.xnnz);
        end
    end
    runhistAS.iter_SNIPAL(1) =  info.iter;
    runhistAS.iter_SSN(1) = info.numSSN;
    runhistAS.Res(1) = Res;
    runhistAS.nnz(1) = info.xnnz;
    subhistAS.n_eachlam(1) = par.n;
    
    ttime = etime(clock,tstart);
    infoAS.iter = 0;
    infoAS.primfeas = info.priminfeas_final;
    infoAS.dualfeas = info.dualinfeas_final;
    infoAS.res_gap = info.res_gap_final;
    infoAS.sigma = info.sigma;
    infoAS.tau = info.tau;
    infoAS.res_kkt = info.res_kkt_final;
    infoAS.ttime = ttime;
    infoAS.totle_iter_SNIPAL = runhistAS.iter_SNIPAL;
    infoAS.totle_iter_SSN = runhistAS.iter_SSN;
    infoAS.Res = Res;
    infoAS.xnnz = info.xnnz;
    
    %% -------------------- Begin the inner LOOP of AS -------------
    parsub = par;
    for itersub = 1:maxitersub
        %% terminal condition
        if Res < tol
            msg = 'AS converges';
            infoAS.termcode = 1;
            x = x_new;
            u = u_new;
            z = z_new;
            break;
        end
        
        %% update x_part and u_new
        ProxATu = Proj_inf(ATu,parsub.lambda);
        cons1 = tol/sqrt(2*lengthI);
        indexJ = (indexI & (abs(ATu-ProxATu) > cons1));
        indexIb = (indexIb | indexJ);
        indexI =  (~indexIb);  
        n_new = sum(indexIb);
        Apart = Ainput(:,indexIb); 
        parsub.n = n_new;
       
        if (iter == 1) && (itersub == 1)
            fields = {'x0','z0','u0'};
            parsub = rmfield(parsub,fields);
        else
            parsub.x0 = x_new(indexIb); parsub.z0 = z_new; parsub.u0 = u_new;
        end
        
        [obj,x_part,z_old,u_new,~,info] = NALM(Apart,b,parsub);
        
        z_tmp = u_new + z_old;
        z_new = z_tmp - Proj_dknorm(z_tmp,parsub.kk,ones(m,1));
        ATu_part = info.ATu;
        x_tmp = x_part-ATu_part;
        x_partnew = x_tmp-Proj_inf(x_tmp,parsub.lambda);
        x_new = zeros(n,1);
        x_new(indexIb) = x_partnew;
        ATu = ATmap(u_new); normATu = norm(ATu);
        Res = norm(ATu + Proj_inf(x_new-ATu,parsub.lambda))/(1+norm(x_new)+normATu);
        
        lengthI = sum(indexI); 
        x = x_new;
        u = u_new;
        z = z_new;
        
        %% Print results of PALM
        if (printyes)
            fprintf('\n %5.1d| %3.2e |%3.2e    %3.2e   %- 3.2e| %- 5.4e %- 5.4e |',...
                itersub,Res,info.priminfeas_final,info.dualinfeas_final,info.res_gap_final,...
                obj(1),obj(2));
            fprintf('%5.1f | %3.2e %3.2e | %3.2e ',info.totletime, info.sigma, info.tau,info.res_kkt_final);
            if Test_stepALM == 1
                fprintf('| %5.1d  %5.1d  %5.1d  %5.1d',info.iter, info.numSSN, parsub.n, info.xnnz);
            end
        end
        runhistAS.iter_SNIPAL(itersub+1) =  info.iter;
        runhistAS.iter_SSN(itersub+1) = info.numSSN;
        runhistAS.Res(itersub+1) = Res;
        runhistAS.nnz(itersub+1) = info.xnnz;
        subhistAS.n_eachlam(itersub+1) =  parsub.n;
    end
    %%
    %% -------------------- END the inner LOOP of AS ----------------------
    %%
    if (itersub == maxitersub)
        msg = ' maximum iteration reached';
        infoAS.termcode = 3;
        x = x_new;
        u = u_new;
        z = z_new;
    end
    
    %----------------------- relative KKT residual ------------------------
    Ax = Amap(x);normu = norm(u);
    eta_x = Res;
    eta_z = norm(u-Proj_dknorm(z+u,par.kk,ones(m,1)))/(1+norm(z)+normu);
    eta_u = norm(Ax-z-b)/(1+normb);
    res_kkt = max([eta_x,eta_z,eta_u]);
    
    %------------------ primal and dual infeasibilities -------------------
    primfeas = eta_u;
    if itersub == 1
        ATu = ATmap(u_new); normATu = norm(ATu);
    end
    ProxATu = Proj_inf(ATu,par.lambda);
    dualfeas1 = norm(ATu-ProxATu)/(1+normATu);
    dualfeas2 = norm(u-Proj_dknorm(u,par.kk,ones(m,1)))/(1+normu);
    dualfeas = max(dualfeas1,dualfeas2);
    
    Axb_tmp = sort(abs(Ax-b),'descend'); normx1 = norm(x,1);
    primobj = sum(Axb_tmp(1:par.kk))+par.lambda*normx1;
    dualobj = -u'*b;
    res_gap = abs(primobj-dualobj)/(1+abs(primobj)+abs(dualobj));
    
    %----------------------the number of nonzeros in x---------------------
    sortx = sort(abs(x),'descend');
    tol_nnzx = 0.999*normx1;
    tmpidex = find(cumsum(sortx) > tol_nnzx);
    if isempty(tmpidex) 
        nnzeros_x = 0;
        eps0 = 0;
    else
        nnzeros_x = tmpidex(1);
        eps0 = sortx(nnzeros_x);
    end
    %----------------------------------------------------------------------
    ttime = etime(clock,tstart);
    Time(iter) = ttime;
    if infoAS.termcode == 1
        infoAS.iter = itersub;
    else
        infoAS.iter = itersub + 1;
    end
    infoAS.primfeas = primfeas;
    infoAS.dualfeas = dualfeas;
    infoAS.res_gap = res_gap;
    infoAS.sigma = info.sigma;
    infoAS.tau = info.tau;
    infoAS.res_kkt = res_kkt;
    infoAS.primobj = primobj;
    infoAS.dualobj = dualobj;
    infoAS.ttime = ttime;
    infoAS.totle_iter_SNIPAL = sum(runhistAS.iter_SNIPAL);
    infoAS.totle_iter_SSN = sum(runhistAS.iter_SSN);
    infoAS.Res = Res;
    infoAS.xnnz = nnzeros_x;
    
    runhistAS.ttime_path(iter) = ttime;
    runhistAS.xnnz_path(iter) = nnzeros_x;
    runhistAS.n_mean(iter) = mean(subhistAS.n_eachlam);
    subhistAS.n_eachlam = [];
    runhistAS.relkkt(iter) = res_kkt;
    runhistAS.relgap(iter) = res_gap;
    runhistAS.primfeas(iter) = infoAS.primfeas;
    runhistAS.dualfeas(iter) = infoAS.dualfeas;
    runhistAS.iterAS(iter) = infoAS.iter;
    runhistAS.iterPAL(iter) = infoAS.totle_iter_SNIPAL;
    runhistAS.iterSSN(iter) = infoAS.totle_iter_SSN;
    runhistAS.termcode(iter) = infoAS.termcode;
    
    X(:,iter) = x;
    U(:,iter) = u;
    Obj(iter,:) = [primobj,dualobj];
    
    if iter == 1
        runhistAS.time(iter) = ttime;
    else
        runhistAS.time(iter) = Time(iter)-Time(iter-1);
    end
    %%
    %% Print Results
    %%
    if (printyes)
        if ~isempty(msg); fprintf('\n %s',msg); end
        fprintf('\n--------------------------------------------------------------');
        fprintf('------------------');
        fprintf('\n  time = %3.2f',ttime);
        fprintf('\n  number iter of AS = %2.0d',itersub);
        fprintf('\n  number iter of SNIPAL = %2.0d',infoAS.totle_iter_SNIPAL);
        fprintf('\n  number iter of SSN = %2.0d',infoAS.totle_iter_SSN);
        fprintf('\n  primobj = %9.8e, dualobj = %9.8e, relgap = %3.2e',primobj,dualobj,infoAS.res_gap);
        fprintf('\n  priminfeas    = %3.2e, dualinfeas    = %3.2e',...
            infoAS.primfeas, infoAS.dualfeas);
        fprintf('\n  relative KKT residual = %3.2e',infoAS.res_kkt);
        fprintf('\n  relative residual of R(x,u) = %3.2e',infoAS.Res);
        fprintf('\n  number of nonzeros in x = %3.0f', infoAS.xnnz);
        fprintf('\n--------------------------------------------------------------');
        fprintf('------------------\n');
    end
    indexI = (abs(x) < eps0);
    indexIb = (abs(x) >= eps0);
    n_new = sum(indexIb);
    Apart = Ainput(:,indexIb);
    
    par.n = n_new;
    par.x0 = x(indexIb); par.z0 = z; par.u0 = u;
    
end
%%%%%%%%%%%----------End Adaptive sieving strategy----------%%%%%%%%%%%%%%%

end

