%%*************************************************************************
%% ADMM_knorm:
%% an alternating direction method of multipliers (ADMM) for solving
%% the following problem
%%     max_{u,v,w}      -<u,b>                                   (D)
%%           s.t.      A^T*u + v = 0,
%%                     w - u = 0, v in B^lam_inf, w in B^1_(k)^*,
%%  where v in R^n, u,w in R^m,
%%  B^lam_inf is a ball with center 0 and radius lambda in the sense 
%%  of ell_inf norm, B^1_(k)^* is a ball with center 0 and radius 1 
%%  in the sense of the duall norm of k-norm.
%% Its dual problem is
%%     min_{x,z}     ||z||_{(k)} + lambda*||x||_1                (P)
%%          s.t.      Ax - z = b,
%% where x in R^n, z in R^m, H is a given positive scalar multiple 
%% of the identity matrix.
%%
%% [obj,x,u,runhist,info] = ADMM_knorm(Ainput,b,options)
%%
%% Input:
%% Ainput, b;
%% options.m = the sample size;
%% options.n = the feature size;
%% options.kk = parameter k in (P);
%% options.lambda = parameter lambda in (P);
%% options.tol = the accuracy tolerance for solving the problem;
%% options.flag_tol = 0, adopt eta_res < tol as the stopping criterion;                      
%%                  = 1, adopt relobj < tol as the stopping criterion;
%% options.obj_opt = the optimal value;
%% options.maxiter = the maximum number of iterations;
%% options.maxtime = the maximum time;
%% options.gamma = the penalty parameter gamma in ADMM;
%% options.tau_ADMM = the step length in ADMM;
%% options.linsolver = 1 : solve linear equation by Cholesky factorization;
%%                     2 : solve linear equation by psqmr;
%%                     3: solve linear equation by semi-proximal dADMM.
%% Output:
%% obj = [Primal objective value, Dual objective value];
%% x = the output primal solution x;
%% u = the output dual solution u;
%% runhist = a structure containing the history of the run;
%% info.iter = the total number of iterations;
%% info.time = total running time;
%% info.time_cpu = total CPU time;
%% info.res_kkt = relative KKT residual;
%% info.eta_res = relative residual based on infeasibilities 
%%                and duality gap.
%% info.relgap = relative duality gap;
%% info.xnnz = the number of nonzero entries for x
%% N-ALM:
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the Appendix D in the supplementary
%% materials of the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%*************************************************************************
function [obj,x,u,runhist,info] = ADMM_knorm(Ainput,b,options)
%%
%% parameter setting
%%
kk = 1;            % the value of parameter k
lambda = 0.1;      % the value of parameter lambda
maxiter = 30000;   % maximum number of iteration
maxtime = 7200;    % the maximum time
tau  = 1.618;      % step length of ADMM
tol = 1e-3;        % the tolerance of ADMM
printyes = 1;      % print information for each iteration of ADMM
printminoryes = 1; % print the initial and final information for ADMM
printprimdualobj_figure = 0; % 1, draw the figure of primal and dual objective values
linsolver = 2;
gamma = 1;         % initial value of panelty parameter gamma
testgamma = 0;     % 1, print and save the information for updating gamma;
                   % otherwise, not
fixlinsolver = 0;  % 0, auto-adjust linsolver according to the dimension of Ainput;
                   % otherwise, fix linsolver
flag_tol = 0; % 1, adopt relobj < tol as the stopping criterion of ADMM;
              % 0, adopt eta_res < tol as the stopping criterion of ADMM.
%-------------------------------------------------------------------------
if isfield(options,'m'); m = options.m; end
if isfield(options,'n'); n = options.n; end
if isfield(options,'kk'); kk = options.kk; end
if isfield(options,'lambda'); lambda = options.lambda; end
if isfield(options,'tol'); tol = options.tol; end
if isfield(options,'flag_tol'), flag_tol = options.flag_tol; end
if flag_tol == 1, obj_opt = options.obj_opt; end
if isfield(options,'maxiter');  maxiter = options.maxiter; end
if isfield(options,'maxtime');  maxtime = options.maxtime; end
if isfield(options,'gamma'); gamma = options.gamma; end
if isfield(options,'tau_ADMM'); tau = options.tau_ADMM; end
if isfield(options,'linsolver'); linsolver = options.linsolver; end
if isfield(options,'printyes'); printyes = options.printyes; end
if isfield(options,'printminoryes'); printminoryes = options.printminoryes; end

% the parameter setting for updating gamma
gamma_mul = 8.9;
gamma_iter = 44;
gamscale1 = 1.8;
gamscale2 = 2.5;
gamscale3 = 1.7;
gamma_siter = 50;
gamma_giter = 500;
if isfield(options,'gamscale1');  gamscale1 = options.gamscale1; end
if isfield(options,'gamscale2');  gamscale2 = options.gamscale2; end
if isfield(options,'gamscale3');  gamscale3 = options.gamscale3; end
if isfield(options,'gamma_mul');  gamma_mul = options.gamma_mul; end
if isfield(options,'gamma_iter');  gamma_iter = options.gamma_iter; end
if isfield(options,'gamma_siter');  gamma_siter = options.gamma_siter; end
if isfield(options,'gamma_giter');  gamma_giter = options.gamma_giter; end
%%
%% Amap and ATmap
%%
tstart = clock;
tstart_cpu = cputime;

A0 = Ainput;
Amap0 = @(x) mexAx(A0,x,0); %Amap0= @(x) A0*x;
ATmap0 = @(y) (y'*A0)'; %AATmap0 = @(z) A0*ATmap0(z);
%%
%% initiallization
%%
normb = norm(b);

u = sparse(m,1); w = u; v = sparse(n,1);
x = v; z = u; u_tmp = zeros(m,1);
%%
%% Scaling
%%
bscale = (max(abs(A0)))';
invbscale = 1./bscale;

Amap = @(x) Amap0(invbscale.*x); % Amap = AB^{-1}
ATmap = @(y) invbscale.*ATmap0(y);
AATmap = @(z) Amap(ATmap(z));
ATAmap = @(w) ATmap(Amap(w));
v = invbscale.*v;
x = bscale.*x;

invB = sparse(1:n,1:n,invbscale);
normA = norm(A0*invB,'fro');
%%
%% linsolver
%%
linnumber = 50000;
if fixlinsolver == 0
    if m < n
        if m <  linnumber
            linsolver = 1;
        else
            linsolver = 2;
        end
    else
        if (n <  linnumber)
            linsolver = 1;
        else
            linsolver = 2;
        end
    end
end

switch(linsolver)
    case 1
        AinvB = A0*invB;
        if m < n
            ImAAT = speye(m) + AinvB*AinvB';
            L = CholHess(ImAAT);
        else
            InATA = speye(n)+AinvB'*AinvB;
            L = CholHess(InATA);
        end
    case 2
        if m < n
            if norm(u) == 0
                ImAATu = sparse(m,1);
            else
                ImAATu = u+AATmap(u);
            end
        else
            InATAu_tmp = sparse(n,1);
        end
    otherwise % linsolver = 3
        AinvB = A0*invB;
        AAT = AinvB*(AinvB)';
        num_eig = 10;
        [P,V] = eigs(AAT,num_eig);
        lam_small = V(end,end);
        V2 = V(1:num_eig-1,1:num_eig-1);
        I2 = speye(num_eig-1);
        Lambda1 = V2-lam_small*I2;
        P2 = P(:,1:end-1);
        Tmap = @(u) lam_small*u+P2*Lambda1*(u'*P2)'-AATmap(u);
        
        v2i2 = diag(V2+I2);
        Lambda2 = 1/(1+lam_small)*I2 - sparse(1:num_eig-1,1:num_eig-1,1./v2i2);
        invTmap = @(rhs) (1/(1+lam_small))*rhs - P2*Lambda2*(rhs'*P2)';
end
%%
%% initial primal and dual infeasibilities
%%
BATu = bscale.*ATmap(u); Av = Amap(v); Bv = bscale.*v;
Rp = Amap(x)-b-z;
Rd1 = BATu+Bv; Rd2 = w - u;

primfeas = norm(Rp)/(1+normb);
dualfeas_1 = norm(Rd1)/(1+max(norm(BATu),norm(Bv)));
dualfeas_2 = norm(Rd2)/(1+max(norm(w),norm(u)));
dualfeas = max(dualfeas_1,dualfeas_2);
maxfeas = max(primfeas,dualfeas);

runhist.cputime(1) = etime(clock,tstart);

%%
%% print initial information
%%
if printyes
    fprintf('\n *******************************************************');
    fprintf('******************************************');
    fprintf('\n \t\t   ADMM  for solving the convex CVaR-based sparse linear regression with k = %4.0d ', kk);
    fprintf('\n ******************************************************');
    fprintf('*******************************************\n');
    if printminoryes
        fprintf('\n problem size: n = %3.0f, nb = %3.0f',n, m);
        fprintf('\n lambda = %g, linsolver = %2.0d, normbarA = %6.3f,  tau = %6.3f',...
            lambda, linsolver, normA,tau);
        fprintf('\n initial primfeas = %3.2e, dualfeas = %3.2e', primfeas, dualfeas);
        fprintf('\n ---------------------------------------------------');
    end
    fprintf('\n  iter|  pinforg    dinforg     relgap |      pobj       dobj    |');
    fprintf(' time |  gamma  ');
    if testgamma == 1
        fprintf('|  tau  [ primfeas   dualfeas   feasratio ]');
    end
end
%%
%% main loop
%%
breakyes = 0;
msg = []; % output the final state of ADMM
%%
for iter = 1:maxiter
    uold = u; wold = w; %vold = v;
    xold = x; zold = z;
    
    %% compute u
    rhsu = (1/gamma)*Rp-Av+wold;
    switch(linsolver)
        case 1     % linsolver = 1: Cholesky factorization
            if m <= n
                u = linsysolvefun(L,rhsu);
            else
                u_tmp = linsysolvefun(L,ATmap(rhsu));
                u = rhsu-Amap(u_tmp);
            end
        case 2     % linsolver = 2: psqmr
            paru.tol = 3*max(0.9*tol,min(1/iter^1.1,0.9*maxfeas));
            if m <= n
                paru.tol = paru.tol*normA;
                paru.tol = 3*normA*max(0.9*tol,min(1/iter^1.1,0.9*maxfeas));
                paru.u0 = uold; paru.printlevel = 0;
                [u,ImAATu,resnrmu,solve_oku] = psqmr_ADMM_new('Matvecmu',AATmap,rhsu,paru,ImAATu);
            else
                paru.u0 = u_tmp; paru.printlevel = 0;
                [u_tmp,InATAu_tmp,resnrmu,solve_oku] = ...
                    psqmr_ADMM_new('Matvecnx',ATAmap,ATmap(rhsu),paru,InATAu_tmp);
                u = rhsu-Amap(u_tmp);
            end
        case 3     % linsolver = 3: semi-proximal ADMM
            rhsu4 = rhsu + Tmap(u);
            u = invTmap(rhsu4);
    end
    ATu = ATmap(u);
    
    %% compute v and w
    v = Proj_inf(xold/gamma-ATu,lambda,invbscale);
    w = Proj_dknorm(u+zold/gamma,kk,ones(m,1));
    Av = Amap(v);
    
    %% update mutilplier x and z
    x = xold - tau*gamma*(ATu+v);
    z = zold - tau*gamma*(w-u);
    Ax = Amap(x);
    
    %% recover primal and dual infesibility  
    Rp = Ax-b-z; BATu = bscale.*ATu; Bv = bscale.*v;
    normu = norm(u); normBATu = norm(BATu);
    primfeas = norm(Rp)/(1+normb);
    Rd1 = BATu+Bv; Rd2 = w-u;
    dualfeas_1 = norm(Rd1)/(1+max(normBATu,norm(Bv)));
    dualfeas_2 = norm(Rd2)/(1+max(norm(w),normu));
    dualfeas = max(dualfeas_1,dualfeas_2);
    maxfeas = max(primfeas,dualfeas);
    
    dualfeas_scale1 = norm(ATu + v)/(1+max(norm(ATu),norm(v)));
    dualfeas_scale = max(dualfeas_scale1,dualfeas_2);
    
    %% record history
    runhist.dualfeas(iter) = dualfeas;
    runhist.primfeas(iter) = primfeas;
    runhist.gamma(iter) = gamma;
    
    if linsolver == 2
        periter_CG = length(resnrmu);
        runhist.psqmruiter(iter) = periter_CG;
    end
    %%---------------------------------------------------------
    %% check for termination
    %%---------------------------------------------------------
    ttime = etime(clock,tstart);
    
    if (maxfeas < 1e-1) || (ttime > maxtime)
        absAxb = sort(abs(Ax-b),'descend');
        primobj = sum(absAxb(1:kk))+lambda*norm(invbscale.*x,1);
        dualobj = -b'*u;
        relgap = abs(primobj-dualobj)/(1+abs(primobj)+abs(dualobj));
        if flag_tol == 0
            eta_res = max(maxfeas,relgap);
            if eta_res < tol
                breakyes = 1;
                msg = 'converged';
            end
        else
            relobj = abs(primobj - obj_opt)/(1 + abs(obj_opt));
            if relobj < tol
                breakyes = 1;
                msg = 'converged';
            end
        end
    end
   
    %%--------------------------------------------------------
    %% print results
    %%--------------------------------------------------------
    if (iter <= 200)
        print_iter = 20;
    elseif (iter <= 2000)
        print_iter = 100;
    else
        print_iter = 200;
    end
    if (rem(iter,print_iter)==1 || iter==maxiter) || (breakyes)
        if (maxfeas >= 1e-1)
            absAxb = sort(abs(Ax-b),'descend');
            primobj = sum(absAxb(1:kk))+lambda*norm(invbscale.*x,1);
            dualobj = -b'*u;
            relgap = abs(primobj-dualobj)/(1+abs(primobj)+abs(dualobj));
        end
        if testgamma == 1
            feasratio = dualfeas_scale/primfeas;
        end
        if (printyes)
            fprintf('\n %5.0d| %3.2e    %3.2e  %- 3.2e| %- 5.4e %- 5.4e |',...
                iter,primfeas,dualfeas,relgap,primobj,dualobj);
            fprintf(' %5.1f| %3.2e',ttime, gamma);
            if testgamma == 1
                fprintf('| %2.3f [ %3.2e    %3.2e   %3.2e ]',tau, primfeas,dualfeas_scale,feasratio);
            end
                       
            if linsolver == 2
                fprintf('[%3.0d %3.2e] %3.0d', periter_CG, paru.tol, solve_oku);
            end
        end
        runhist.primobj(iter)   = primobj;
        runhist.dualobj(iter)   = dualobj;
        runhist.time(iter)      = ttime;
        runhist.relgap(iter)    = relgap;
    end
    ttime = etime(clock,tstart);
    if (breakyes > 0) || (ttime > maxtime)
        fprintf('\n   breakyes = %3.1f,',breakyes);
        if flag_tol == 1
            fprintf(' relobj = %3.2e',relobj);
        else
            fprintf(' eta_res = %3.2e',eta_res);
        end
        break;
    end
    %%-----------------------------------------------------------
    %% update penalty parameter gamma
    feasratio = dualfeas_scale/primfeas;
    maxfeasratio = max(feasratio,1/feasratio);
    
    if testgamma == 1
        runhist.feasratio(iter) = feasratio;
        runhist.maxfeasratio(iter) = maxfeasratio;
    end
    
    if maxfeasratio <= gamma_siter
        gammascale = gamscale1;
    elseif maxfeasratio > gamma_giter
        gammascale = gamscale2;
    else
        gammascale = gamscale3;
    end
    
    if (rem(iter,gamma_iter) == 0)
        if feasratio > gamma_mul
            gamma = gamma*gammascale;
        elseif feasratio < 1/gamma_mul
            gamma = gamma/gammascale;
        end
    end
    
    
end
%------------------------------End dADMM main loop-------------------------
%%-----------------------------------------------------------------
%% recover original variables
%%-----------------------------------------------------------------
if (iter == maxiter)
    msg = ' maximum iteration reached';
    info.termcode = 3;
end

x = invbscale.*x;
%v = bscale.*v;
Axb = Amap0(x)-b;
ATu = ATmap0(u);

absAxb = sort(abs(Axb),'descend');
primobj = lambda*norm(x,1)+sum(absAxb(1:kk));
dualobj = -u'*b;
obj = [primobj, dualobj];
relgap = abs(primobj-dualobj)/(1+abs(primobj)+abs(dualobj));
eta_res = max(maxfeas,relgap);

etax = norm(ATu+Proj_inf(x-ATu,lambda))/(1+norm(x)+norm(ATu));
etaz = norm(u-Proj_dknorm(u+z,kk,ones(m,1)))/(1+normu+norm(z));
etau = primfeas;

eta_kkt = max([etax,etaz,etau]);
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

runhist.m = m;
runhist.n = n; 
ttime = etime(clock,tstart);
ttime_cpu = cputime - tstart_cpu;
runhist.iter = iter;
runhist.totaltime = ttime;
runhist.primobjorg = primobj;
runhist.dualobjorg = dualobj;
runhist.maxfeas = maxfeas;

if linsolver == 2
    ttCG = sum(runhist.psqmruiter); % the total iterations of CG
    info.ttCG = ttCG;
end

info.relgap = relgap;
info.iter = iter;
info.time = ttime;
info.time_cpu = ttime_cpu;
info.gamma = gamma;
info.res_kkt = eta_kkt;
info.xnnz = nnzeros_x;
info.eta_res = eta_res;
if flag_tol == 1
    info.relobj = relobj;
end

if (printminoryes)
    if ~isempty(msg); fprintf('\n   %s',msg); end
    fprintf('\n--------------------------------------------------------------');
    fprintf('------------------');
    fprintf('\n  number iter = %2.0d',iter);
    fprintf('\n  time = %3.2f',ttime);
    fprintf('\n  time per iter = %5.4f',ttime/iter);
    fprintf('\n  cputime = %3.2f', ttime_cpu);
    fprintf('\n  primobj = %9.8e, dualobj = %9.8e, relgap = %3.2e',primobj,dualobj, relgap);
    fprintf('\n  primfeas    = %3.2e, dualfeas    = %3.2e, relkkt  = %3.2e',...
        primfeas, dualfeas, eta_kkt);
    fprintf('\n  eta_res  = %3.2e',eta_res);
    fprintf('\n  number of nonzeros in x (abs(x)>1e-10) = %3.0f', nnzeros_x);
    
    if linsolver == 2
        fprintf('\n  Total CG number = %3.0d, CG per iter = %3.1f', ttCG, ttCG/iter); % notice
    end
    
    if flag_tol == 1
        fprintf('\n  relobj = %3.2e',relobj);
    end
    fprintf('\n--------------------------------------------------------------');
    fprintf('------------------\n');
end

if printprimdualobj_figure == 1
    figure
    iter_sub = find(runhist.primobj);
    iter_len = length(iter_sub);
    iter_t = iter_sub(10:iter_len);
    plot(iter_t,runhist.primobj(iter_t),'--*r');
    hold on
    plot(iter_t,runhist.dualobj(iter_t),'-.pb');
    legend('primal objective value','dual objective value');
    xlabel('iter');
    title('The primal and dual objective values');
end

runhist.primobj(iter)   = primobj;
runhist.dualobj(iter)   = dualobj;




