%%*************************************************************************
%% NALM:
%% A semismooth-Newton based proximal augmented Lagrangian method for solving
%% the convex CVaR-based sparse linear regression:
%%
%% (P)   minimize_{x in R^n} {||A*x - b||_(k) + lambda ||x||_1}
%%
%% (D)   maximize_{u in R^m}         -<u,b>
%%       subject to         ||A^T*u||_{infty} <= lambda,
%%                                ||u||_(k)^* <= 1}.
%% where A in R^{m * n}, b in R^m and lambda > 0 are given.
%%
%% [obj,x,u,z,runhist,info] = N_ALM(Ainput,b,OPTIONS)
%% Input:
%% Ainput,b;
%% OPTIONS.tol = the accuracy tolerance for solving the problem;
%% OPTIONS.maxiter = the maximum number of outer iteration;
%% OPTIONS.maxitersub = the maximum number of inner iteration;
%% OPTIONS.maxtime = the maximum time of the N-ALM;
%% OPTIONS.lambda = the value of parameter lambda in (P);
%% OPTIONS.kk = the value of parameter k in (P);
%% OPTIONS.sigma = the initial value of sigma in the proximal ALM;
%% OPTIONS.tau = the initial value of tau in the proximal ALM;
%% OPTIONS.sigmascale = the positive number for updating sigma;
%% OPTIONS.tauscale = the positive number for updating tau;
%% OPTIONS.flag_tol = 0, adopt eta_res < tol as the stopping criterion of the N-ALM.
%%                  = 1, adopt  relobj < tol as the stopping criterion of the N-ALM;
%%                  = 2, adopt  relkkt < tol as the stopping criterion of the N-ALM.
%% Output:
%% obj = [Primal objective value, Dual objective value];
%% x = the primal variable;
%% z = the primal slack variable;
%% u = the dual variable;
%% runhist = a structure containing the history of the run;
%% info.iter = the total number of outer iterations;
%% info.numSSN = the total number of inner iterations;
%% info.priminfeas_final = primal infeasibility;
%% info.dualinfeas_final = dual infeasibility;
%% info.res_kkt_final = relative KKT residual;
%% info.eta_res = relative residual based on infeasibilities and duality gap.
%% info.totletime = total time;
%% info.xnnz = the number of nonzero entries for x
%% N-ALM:
%% Copyright (c) 2022 by
%% Can Wu, Ying Cui, Donghui Li, Defeng Sun
%%*************************************************************************

function [obj,x,z,u,runhist,info] = NALM(Ainput,b,OPTIONS)
%%
%% Input parameters
%%
tol = 1e-6;
kk = 1;
lamc = 1;
sigma = 0.01;
tau = 1;
printyes = 1; % print in ALM
printlevel = 1; % print in SSN
flag_tol = 2;
maxiter = 200; % the maximum iteration numbers of the ALM
maxitersub = 200; % the maximum iteration numbers of the SSN
maxtime = 7200; % the maximum time of the N-ALM

if isfield(OPTIONS,'tol'), tol = OPTIONS.tol; end
if isfield(OPTIONS,'maxiter'), maxiter = OPTIONS.maxiter; end
if isfield(OPTIONS,'maxitersub'), maxitersub = OPTIONS.maxitersub; end
if isfield(OPTIONS,'maxtime'), maxtime = OPTIONS.maxtime; end
if isfield(OPTIONS,'tau'), tau = OPTIONS.tau; end
if isfield(OPTIONS,'sigma'), sigma = OPTIONS.sigma; end
if isfield(OPTIONS,'sigmascale'), sigmascale = OPTIONS.sigmascale; end
if isfield(OPTIONS,'tauscale'), tauscale = OPTIONS.tauscale; end
if isfield(OPTIONS,'kk'), kk = OPTIONS.kk; end
if isfield(OPTIONS,'m'), m = OPTIONS.m; end
if isfield(OPTIONS,'n'), n = OPTIONS.n; end
if isfield(OPTIONS,'lambda'), lambda = OPTIONS.lambda;else, lambda = lamc*max(abs(b'*Ainput)); end
if isfield(OPTIONS,'x0'), x0 = OPTIONS.x0; else, x0 = zeros(n,1); end
if isfield(OPTIONS,'z0'), z0 = OPTIONS.z0; else, z0 = zeros(m,1); end
if isfield(OPTIONS,'u0'), u0 = OPTIONS.u0; else, u0 = zeros(m,1); end
if isfield(OPTIONS,'printyes'), printyes = OPTIONS.printyes; end
if isfield(OPTIONS,'printlevel'), printlevel = OPTIONS.printlevel; end
if isfield(OPTIONS,'flag_tol'), flag_tol = OPTIONS.flag_tol; end
if flag_tol == 1, obj_opt = OPTIONS.obj_opt; end

par.kk = kk;
par.m = m;
par.n = n;
par.lambda = lambda;
par.sigma = sigma;
par.tau = tau;
par.indexJ = false(n,1);
par.AJ = [];
par.AJTAJ = [];
par.num_T = 0;

if printlevel == 1
    printsub = 1; % print the information of stagnation in SSN
else
    printsub = 0;
end

breakyes = 0;
numSSN = 0;
numCG = 0; % run CG
iterCG = 0; % the totle iterations of CG
msg = [];

%testgrad = 0;
%testJacobian = 0;

direc = zeros(m,1);
Im = ones(m,1);
par.Im = Im;
tiny = 1e-10;
%%
%% Amap and ATmap
%%
if isstruct(Ainput)
    if isfield(Ainput,'A'); A = Ainput.A; end
    if isfield(Ainput,'Amap'); Amap = Ainput.Amap; end
    if isfield(Ainput,'ATmap'); ATmap = Ainput.ATmap; end
else
    A = Ainput;
    Amap = @(x) mexAx(A,x,0);
    ATmap = @(y) (y'*A)'; %ATmap = @(x) mexAx(A,x,1);
end

%%
%% Generate initial point
%%
tstart = clock;
tstart_cpu = cputime;

if printyes
    fprintf('\n *******************************************************');
    fprintf('******************************************');
    fprintf('\n \t\t   The N-ALM  for solving the convex CVaR-based sparse linear regression with k = %6.3f', par.kk);
    fprintf('\n ******************************************************');
    fprintf('*******************************************\n');
    fprintf('\n problem size: m = %3.0f, n = %3.0f,',m, n);
    fprintf(' lambda = %g, sigma0 = %g,  tau0 = %g', par.lambda, par.sigma, par.tau);
    fprintf('\n sigmascale = %3.2f, tauscale = %3.2f', sigmascale, tauscale);
    fprintf('\n ---------------------------------------------------');
    fprintf('---------------------------------------------');
    fprintf('\n  iter|  priminfeas  dualinfeas   relgap |    pobj       dobj      |');
    fprintf(' time |  sigma    tau     tau/sig  | res_kkt  |');
end
%%
%% print the initial information
%%
iter = 0; x_new = x0; u_new = u0; z_new = z0;
normx = norm(x_new); normu = norm(u_new);
normz = norm(z_new); normb = norm(b);
if normx == 0
    Ax = sparse(m,1);
else
    Ax = Amap(x_new);
end
if normu == 0
    ATu = sparse(n,1);
else
    ATu = ATmap(u_new);
end
normATu = norm(ATu);

eta_x = norm(ATu+Proj_inf(x_new-ATu,par.lambda))/(1+normx+normATu);
eta_z = norm(u_new-Proj_dknorm(z_new+u_new,par.kk,Im))/(1+normz+normu);
eta_u = norm(Ax-z_new-b)/(1+normb);
res_kkt = max([eta_x,eta_z,eta_u]);

primfeas = eta_u;
dualfeas1 = norm(ATu-Proj_inf(ATu,par.lambda))/(1+normATu);
dualfeas2 = norm(u_new-Proj_dknorm(u_new,par.kk,Im))/(1+normu);
dualfeas = max(dualfeas1,dualfeas2);

Axb_tmp = sort(abs(Ax-b),'descend');
primobj = sum(Axb_tmp(1:par.kk))+par.lambda*norm(x_new,1);
dualobj = -u_new'*b;
res_gap = abs(primobj-dualobj)/(1+abs(primobj)+abs(dualobj)); % the duality gap

w1 = x_new - par.sigma*ATu;
xp =  Proj_inf(w1,par.sigma*par.lambda);
w2 = z_new + par.sigma*u_new;
zp = Proj_dknorm(w2,par.kk*par.sigma,par.sigma*Im);

ttime = etime(clock,tstart);
if (printyes)
    fprintf('\n %5.1d|  %3.2e    %3.2e   %- 3.2e| %- 5.4e %- 5.4e |',...
        iter,primfeas,dualfeas,res_gap,primobj,dualobj);
    fprintf('%5.1f | %3.2e %3.2e %3.2e |',ttime, par.sigma, par.tau, par.tau/par.sigma);
    fprintf(' %3.2e |',res_kkt);
end
%%
%% Proximal ALM main
%%
for iter = 1:maxiter
    par.iter = iter;
    %%
    %% SSN --> u_new
    %%
    w1_sub = w1; xp_sub = xp; w2_sub = w2; zp_sub = zp;
    w2_zp_sub = w2_sub-zp_sub; w1_xp_sub = w1_sub-xp_sub;
    gradphi = Amap(w1_xp_sub)-w2_zp_sub-b; % -gradphi
    phi = -(1/(2*par.sigma))*(norm(w1_xp_sub)^2+norm(w2_zp_sub)^2)-u_new'*b; % -phi
    u_snew = u_new; x_snew = x_new; z_snew = z_new;
    parsub = par; normdelu = 0; subhist.solve_ok = [];
    dualfeas_sub = dualfeas; primfeas_sub = primfeas;
    cntAmap_sub = 0; cntATmap_sub = 0;
    ATu_snew = ATu;
    subhist.r_indexJ = [];
    for itersub = 1:maxitersub
        break_ok = 0;
        %% Stopping criterion: (A) and (B)
        normgradphi = norm(gradphi);
        subhist.gradphi(itersub) = normgradphi;
        const1 = min(sqrt(parsub.tau),1)/parsub.sigma;
        epsilon = 1/(iter^(1.2));
        delta = min(1/(iter^(1.2)),0.5);
        tolsub1 = const1*epsilon;
        if itersub == 1
            tolsub2 = 0;
        else
            deluxzsq = sqrt(parsub.tau*normdelu^2+norm(x_snew-x_new)^2+norm(z_snew-z_new)^2);
            tolsub2 = const1*delta*deluxzsq;
        end
        tolsub = min([1,tolsub1,tolsub2]);
        
        if (itersub > 1)
            if (normgradphi <= tolsub)  || (normgradphi < 1e-12)
                if printlevel
                    msg = ' good termination for ALM:';
                    fprintf('\n       %s  ',msg);
                    fprintf('\n        iters=%2.0d,  gradphi=%3.2e, tolsub=%3.2e, const1=%3.2e',...
                        itersub, normgradphi, tolsub, const1);
                end
                break_ok = -1;
                break;
            end
        end
        if (itersub > 50)
            ratio_gradphi = subhist.gradphi(itersub-9:itersub)./subhist.gradphi(itersub-10:itersub-1);
            if (min(ratio_gradphi) > 0.997) && (max(ratio_gradphi) < 1.003)
                break_ok = -2;
                if (printlevel);  fprintf('stagnate'); end
                break;
            end
        end
        %% Choose linsolver
        [flag_case,r_tmp,parsub] = Generatedash_Jacobi(A,w1_sub,w2_sub,parsub);
        upnumber = 10000; r_indexJ = parsub.r_indexJ;
        subhist.r_indexJ(itersub) = r_indexJ;
        if (m <= upnumber)
            if (m <= 1000)
                linsolver = 'd_direct';
            elseif (r_indexJ <= max(0.01*n, upnumber))
                linsolver = 'd_direct';
            end
        end
        if ((r_indexJ <= m) && (r_indexJ <= upnumber))
            linsolver = 'p_direct';
        end
        if (m > 5e3 && r_indexJ >= 3000) || (m > 2000 && r_indexJ > 8000) || (m > 100 && r_indexJ > 1e4)
            linsolver = 'd_pcg';
        end
        
        %% prepare for psqmr
        if strcmp(linsolver,'d_pcg')
            % -------------------choose maxitpsqmr----------------------------
            maxitpsqmr = 100; % adjusting
            if ((dualfeas_sub > 1e-3) || (itersub <= 5))
                maxitpsqmr = max(maxitpsqmr,100);
            elseif (dualfeas_sub > 1e-4)
                maxitpsqmr = max(maxitpsqmr,300);
            elseif (dualfeas_sub > 1e-5)
                maxitpsqmr = max(maxitpsqmr,400);
            elseif (dualfeas_sub > 5e-6)
                maxitpsqmr = max(maxitpsqmr,500);
            end
            parsub.minitpsqmr = 3;
            %-----------------------choose stagnate_check_psqmr------------
            stagnate_check_psqmr = 20; % adjusting
            if (dualfeas_sub > 1e-4)
                stagnate_check_psqmr = max(stagnate_check_psqmr,20);
            else
                stagnate_check_psqmr = max(stagnate_check_psqmr,60);
            end
            if (itersub > 3) && (dualfeas_sub < 5e-5) && (length(subhist.solve_ok)>=3)
                if all(subhist.solve_ok(itersub-(1:3)) <= -1)
                    stagnate_check_psqmr = max(stagnate_check_psqmr,100);
                end
            end
            parsub.stagnate_check_psqmr = stagnate_check_psqmr;
            %------------------------- tolpsqmr -------------------------------
            rhs = gradphi;
            if (itersub > 1)
                prim_ratio = primfeas_sub/subhist.primfeas(itersub-1);
                dual_ratio = dualfeas_sub/subhist.dualfeas(itersub-1);
            else
                prim_ratio = 0; dual_ratio = 0;
            end
            tolpsqmr = min(5e-3, 0.1*norm(rhs));
            const2 = 1;
            if itersub > 1 && (prim_ratio > 0.5 || primfeas > 0.1*subhist.primfeas(1))
                const2 = 0.5*const2;
            end
            if (dual_ratio > 1.1); const2 = 0.5*const2; end
            tolpsqmr = const2*tolpsqmr;
            L_pre.invdiagM = [];
        end
        
        %% compute Newton direction and step length
        rhs = gradphi;
        switch linsolver
            case 'd_direct'
                if r_tmp == 0    % no AJAJT and W
                    direc = (parsub.sigma/parsub.tau)*rhs;
                else
                    W = Cumpute_matrix_W(flag_case,parsub);
                    Hs = parsub.sigma*(parsub.AJ*parsub.AJ' + W) + (parsub.tau/parsub.sigma)*speye(m);
                    if m <= 1500
                        direc = Hs\rhs;
                    else
                        L = CholHess(Hs);
                        direc = linsysolvefun(L,rhs);
                    end
                end
                resnrm = 0; solve_ok = 1;
            case 'p_direct'
                if (parsub.r_indexJ == 0)
                    [direc,~,~] = invDvec(rhs,flag_case,parsub,0);
                else
                    [invDrhs,T,parsub] = invDvec(rhs,flag_case,parsub,1);
                    if parsub.r_indexJ <= 1500
                        d_tmp = T\((invDrhs)'*parsub.AJ)';
                    else
                        L = CholHess(T);
                        d_tmp = linsysolvefun(L,((invDrhs)'*parsub.AJ)');
                    end
                    [invDrhs_tmp,~,~] = invDvec(parsub.AJ*d_tmp,flag_case,parsub,0);
                    direc = invDrhs - invDrhs_tmp;
                    parsub.num_T = parsub.num_T + 1;
                end
                resnrm = 0; solve_ok = 1;
            case 'd_pcg'
                parsub.d = direc;
                C = Cumpute_matrix_C(flag_case,parsub);
                [direc,resnrm,solve_ok] = psqmr_knorm_N_ALM('matvec_N_ALM',C,rhs,parsub,L_pre,...
                    tolpsqmr,maxitpsqmr,0);
                numCG = numCG + 1;
        end
        subhist.solve_ok(itersub) = solve_ok;
        psqmriter = length(resnrm)-1; % nuber of iterations of psqmr        iterCG =[iterCG psqmriter];
        subhist.psqmr(itersub) = psqmriter;
        
        %% Strongly Wolfe search for stepsize: update u_snew
        if (itersub<=3) && (dualfeas_sub > 1e-4) || (itersub <3)
            stepop = 1;
        else
            stepop = 2;
        end
        steptol = 1e-7; step_op.stepop=stepop;
        [u_snew,x_snew,z_snew,alpha,iterstep,g0,w1_sub,w2_sub,phi,maxiterfs,ATdu,normdelu]= findstep...
            (parsub,ATmap,b,w1_sub,w2_sub,x_snew,z_snew,u_snew,u_new,direc,gradphi,phi,steptol,step_op);
        numSSN = numSSN + 1;
        cntATmap_sub = cntATmap_sub + 1;
        if alpha < tiny; break; end %if alpha < tiny; breakyes =11; break; end
        
        %% update primfeas_sub, dualfeas_sub and -gradphi
        ATu_snew = ATu_snew + alpha*ATdu; %ATu_snew = ATmap(u_snew);
        Axb_sub = Amap(x_snew)-b;
        cntAmap_sub = cntAmap_sub +  1;
        
        primfeas_tmp = Axb_sub-z_snew;
        primfeas_sub = norm(primfeas_tmp)/(1+normb);
        normATu_snew = norm(ATu_snew); normu_snew = norm(u_snew);
        dualfeas1_sub = norm(ATu_snew-Proj_inf(ATu_snew,parsub.lambda))/(1+normATu_snew);
        dualfeas2_sub = norm(u_snew-Proj_dknorm(u_snew,parsub.kk,Im))/(1+normu_snew);
        dualfeas_sub = max(dualfeas1_sub,dualfeas2_sub);
        subhist.primfeas(itersub) = primfeas_sub;
        subhist.dualfeas(itersub) = dualfeas_sub;
        
        subhist.solve_ok(itersub) = solve_ok;
        gradphi = primfeas_tmp-(parsub.tau/parsub.sigma)*(u_snew-u_new);
        
        %% check for stagnation
        if (itersub > 4)
            idx = [max(1,itersub-3):itersub];
            tmp = subhist.primfeas(idx);
            ratio = min(tmp)/max(tmp);
            if (all(subhist.solve_ok(idx) <= -1)) && (ratio > 0.9) ...
                    && (min(subhist.psqmr(idx)) == max(subhist.psqmr(idx))) ...
                    && (max(tmp) < 5*tol)
                fprintf('#')
                break_ok = 1;
                break;
            end
            const3 = 0.7;
            priminf_1half  = min(subhist.primfeas(1:ceil(itersub*const3)));
            priminf_2half  = min(subhist.primfeas(ceil(itersub*const3)+1:itersub));
            priminf_best   = min(subhist.primfeas(1:itersub-1));
            priminf_ratio  = subhist.primfeas(itersub)/subhist.primfeas(itersub-1);
            %dualinf_ratio  = subhist.dualfeas(itersub)/subhist.dualfeas(itersub-1);
            stagnate_idx   = find(subhist.solve_ok(1:itersub) <= -1);
            stagnate_count = length(stagnate_idx);
            idx2 = [max(1,itersub-7):itersub];
            if (itersub >= 10) && all(subhist.solve_ok(idx2) == -1) ...
                    && (priminf_best < 1e-2) && (dualfeas_sub < 1e-3)
                tmp = subhist.primfeas(idx2);
                ratio = min(tmp)/max(tmp);
                if (ratio > 0.5)
                    if (printsub); fprintf('##'); end
                    break_ok = 2;
                end
            end
            if (itersub >= 15) && (priminf_1half < min(2e-3,priminf_2half)) ...
                    && (dualfeas_sub < 0.8*subhist.dualfeas(1)) && (dualfeas_sub < 1e-3) ...
                    && (stagnate_count >= 3)
                if (printsub); fprintf('###'); end
                break_ok = 3;
            end
            if (itersub >= 15) && (priminf_ratio < 0.1) ...
                    && (primfeas_sub < 0.8*priminf_1half) ...
                    && (dualfeas_sub < min(1e-3,2*primfeas_sub)) ...
                    && ((primfeas_sub < 2e-3) || (dualfeas_sub < 1e-5 && primfeas_sub < 5e-3)) ...
                    && (stagnate_count >= 3)
                if (printsub); fprintf(' $$'); end
                break_ok = 4;
            end
            
        end
        %% Print results of SSN
        if (itersub <= 10)
            print_itersub = 1;
        else
            print_itersub = 1;
        end
        if strcmp(linsolver,'d_direct')
            solver = 1;
        elseif strcmp(linsolver,'p_direct')
            solver = 2;
        elseif strcmp(linsolver,'d_pcg')
            solver = 3;
        end
        
        runhist.phi(itersub) = -phi;
        if (printlevel) && (rem(itersub,print_itersub)==0 || itersub == maxitersub)
            fprintf('\n\t\t [%2.0d]   [%3.2e  %3.2e] [%3.2e  %3.2e]',...
                itersub, normgradphi, tolsub, primfeas_sub, dualfeas_sub);
            fprintf(' flag=%1.1f, g0=%3.2e, lin=%2.0d', flag_case, -g0, solver);
            fprintf(' [%3.2e %2.0f] ',alpha,iterstep);
            %------------test----------------------------------------------
            fprintf(' [%3.0f %3.0f %3.0f %3.0f] ', parsub.r_indexJ,...
                parsub.index_alpha,parsub.index_beta,parsub.index_gamma);
            %             if itersub > 1
            %                 delphi = runhist.phi(itersub)-runhist.phi(itersub-1);
            %                 fprintf(' delphi=%3.2e ', delphi);
            %             end
            %-----------end------------------------------------------------
            if strcmp(linsolver,'d_pcg') && ~(parsub.r_indexJ==0)
                fprintf('[%3.1e %3.1e %3.0d]',tolpsqmr,resnrm(end),psqmriter);
            end
            if maxiterfs == 1
                fprintf('$');
            end
        end
        if (break_ok > 0)
            break;
        end
        
    end
    %% End SSN
    runhist.r_indexJ(iter)= mean(subhist.r_indexJ);
    par.AJ = parsub.AJ; par.AJTAJ = parsub.AJTAJ; par.indexJ = parsub.indexJ;
    par.num_T = parsub.num_T; ATu = ATu_snew;
    if break_ok  < 0
        runhist.iterSSN(iter) = itersub - 1;
    else
        runhist.iterSSN(iter) = itersub;
    end
    %% update u_new, z_new, x_new
    u_new = u_snew;
    x_new = x_snew;
    z_new = z_snew;
    ATu_new = ATu_snew;
    primfeas = primfeas_sub;
    dualfeas = dualfeas_sub;
    Axb_new = Axb_sub;
    %% Stopping criterion: relative duality gap
    
    Axb_tmp = sort(abs(Axb_new),'descend');
    k_norm_Axb = sum(Axb_tmp(1:par.kk)); laml1norm  = par.lambda*norm(x_new,1);
    primobj =  k_norm_Axb + laml1norm;
    dualobj = -u_new'*b;
    res_gap = abs(primobj-dualobj)/(1+abs(primobj)+abs(dualobj));
    eta_res = max([primfeas,dualfeas,res_gap]);
    
    normATu_new = normATu_snew;
    normu_new = normu_snew;
    eta_x = norm(ATu_new+Proj_inf(x_new-ATu_new,par.lambda))/(1+norm(x_new)+normATu_new);
    eta_z = norm(u_new-Proj_dknorm(z_new+u_new,par.kk,Im))/(1+norm(z_new)+normu_new);
    eta_u = primfeas;
    res_kkt = max([eta_x,eta_z,eta_u]);
    
    info.res_gap(iter) = res_gap;
    info.priminfeas(iter) = primfeas;
    info.dualinfeas(iter) = dualfeas;
    runhist.primobj(iter) = primobj;
    runhist.dualobj(iter) = dualobj;
    runhist.ratio(iter) = dualfeas/primfeas;
    runhist.res_kkt(iter) = res_kkt;
    runhist.eta_x(iter) = eta_x;
    runhist.eta_z(iter) = eta_z;
    runhist.eta_u(iter) = eta_u;
    
    switch(flag_tol)
        case 0
            if eta_res < tol
                breakyes = 1;
                msg = ' converged';
            end
        case 1
            relobj = abs(primobj - obj_opt)/(1 + abs(obj_opt));
            if relobj < tol
                breakyes = 1;
                msg = ' converged';
            end
        case 2
            if res_kkt < tol
                breakyes = 1;
                msg = ' converged';
            end
    end
    
    %% Print results of PALM
    if (iter <= 20)
        print_iter = 1;
    elseif (iter <= 50)
        print_iter = 1;
    else
        print_iter = 1;
    end
    ttime = etime(clock,tstart);
    runhist.ttime(iter) = ttime;
    if (rem(iter,print_iter)==0 || iter==maxiter) || (breakyes==1) || (ttime > maxtime)
        if (printyes)
            fprintf('\n %5.0d|  %3.2e    %3.2e   %- 3.2e| %- 5.4e %- 5.4e |',...
                iter,primfeas,dualfeas,res_gap,primobj,dualobj);
            fprintf('%5.1f | %3.2e %3.2e %3.2e |',ttime, par.sigma, par.tau, par.tau/par.sigma);
            fprintf(' %3.2e |',res_kkt);
        end
        
    end
    
    %% update sigma and tau
    if  mod(iter,3) == 0
        par.sigma = min(par.sigma*sigmascale,1e6);
        par.tau = max(par.tau/tauscale,1e-6);
    end
    
    %%
    subinfo.cntAmap(iter) = cntAmap_sub + 1;
    subinfo.cntATmap(iter) = cntATmap_sub; % the same as numSSN
    %% termination
    if (breakyes > 0) || (ttime > maxtime) || (iter == maxiter)
        if (printyes)
            fprintf('\n  breakyes = %3.1f, %s,',breakyes,msg);
            switch(flag_tol)
                case 0
                    fprintf(' eta_res = %3.2e',eta_res);
                case 1
                    fprintf(' relobj = %3.2e',relobj);
                case 2
                    fprintf(' relkkt = %3.2e',res_kkt);
            end
        end
        x = x_new; z = z_new; u = u_new;
        break;
    end
    
    w1 = x_new - par.sigma*ATu_new;
    xp =  Proj_inf(w1,par.sigma*par.lambda);
    w2 = z_new + par.sigma*u_new;
    zp = Proj_dknorm(w2,par.kk*par.sigma,par.sigma*Im);
end
%% End PALM
if (iter == maxiter)
    msg = ' maximum iteration reached';
    info.termcode = 3;
end
if (ttime > maxtime)
    msg = ' maximum time reached';
    info.termcode = 4;
end
ttime = etime(clock,tstart);
ttime_cpu = cputime - tstart_cpu;
obj = [primobj, dualobj];

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
%--------------------------------------------------------------------------



info.minx = min(x); info.maxx = max(x);
info.totle_cntAmap = sum(subinfo.cntAmap);
info.totle_cntATmap = sum(subinfo.cntATmap);
%%
%% Print Results
%%
if (printyes)
    if ~isempty(msg); fprintf('\n %s',msg); end
    fprintf('\n--------------------------------------------------------------');
    fprintf('------------------');
    fprintf('\n  number iter = %2.0d',iter);
    fprintf('\n  number iter of SSN = %2.0d',numSSN);
    fprintf('\n  SSN per iter = %3.0f', numSSN/iter);
    fprintf('\n  average number of PSQMR = %3.1f',mean(iterCG));
    fprintf('\n  time = %3.2f',ttime);
    fprintf('\n  time per iter = %5.4f',ttime/iter);
    fprintf('\n  cputime = %3.2f', ttime_cpu);
    fprintf('\n  primobj = %9.8e, dualobj = %9.8e, relgap = %3.2e',primobj,dualobj, res_gap);
    fprintf('\n  priminfeas    = %3.2e, dualinfeas    = %3.2e',...
        primfeas, dualfeas);
    fprintf('\n  relative KKT residual = %3.2e',res_kkt);
    fprintf('\n  knorm(Ax-b) = %5.4f, lambda*norm(x,1) = %5.4f',k_norm_Axb,laml1norm);
    fprintf('\n  number of nonzeros in x = %3.0f', nnzeros_x);
    fprintf('\n  min(x) = %3.2e, max(x) = %3.2e', info.minx, info.maxx);
    fprintf('\n  Amap cnt = %4.0f, ATmap cnt = %4.0f', info.totle_cntAmap,info.totle_cntATmap);
    if strcmp(linsolver,'d_pcg')
        fprintf('\n  Total CG number = %3.0d, CG per iter = %3.1f', numCG, numCG/iter); % notice
    end
    fprintf('\n--------------------------------------------------------------');
    fprintf('------------------\n');
end

%%
%% record history
%%
info.x = x;
info.z = z;
info.u = u;
info.iter = iter;
info.numSSN = numSSN;
info.res_kkt_final = res_kkt;
info.res_gap_final = res_gap;
info.priminfeas_final = primfeas;
info.dualinfeas_final = dualfeas;
info.eta_res = eta_res;
info.xnnz = nnzeros_x;
info.totletime = ttime;
info.sigma = par.sigma;
info.tau = par.tau;
info.ATu = ATu_new;
info.z = z;
info.breakyes = breakyes;
if flag_tol == 1
    info.relobj = relobj;
end
info.index_nnz = index_nnz;
info.index_0 = index_0;

end
