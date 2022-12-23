%%*************************************************************************
%% dPPA_SSN:
%% A proximal point algorithm for solving the dual problem of the 
%% following subproblem of MM algorithm
%%
%% minimize    { ||A*x - b||_(k1) + lambda ||x||_1 - ||A*x_old - b||_(k2)  
%% x in R^n       - < a, A*x - A*x_old > + (rho/2)||x - x_old||^2 
%%                + (sigma/2)||A*x - A*x_old||^2 }
%% where a is a subgradient of f_2(x) := ||Ax-b||_(k2) at x = x_old,
%% rho > 0, sigma > 0, x_old is a given vector.
%%
%% [x,u,z,info,runhist] = dPPA_SSN(A,b,OPTIONS)
%% Input: A, b
%% OPTIONS.a = the subgradient of f_2 at x_old;
%% OPTIONS.rho = parameter rho in the above subproblem;
%% OPTIONS.sigma = parameter sigma in the above subproblem;
%% OPTIONS.eta0 = the initial value of the parameter eta in PPA;
%% OPTIONS.etascale = coefficient in [1,+inf) for updating eta;
%% OPTIONS.obj_gap_mm = previous relative residual for two adjacent
%%                      objectives in MM;
%% OPTIONS.x0 = the initial point of primal variable x;
%% OPTIONS.u0 = the initial point of dual variable u;
%% OPTIONS.printPPA = 1, print in PPA; 0, does not print in PPA;
%% OPTIONS.printSSN = 1, print in SSN; 0, does not print in SSN;
%% Output:
%% x = the output primal solution x;
%% u = the output dual solution u;
%% z = the output primal solution z;
%% info.err_ppa = the residual of PPA;
%% info.tol_ppa = (sigma/2)||A*(x - x_old)||^2;
%% info.primfeas = the primal infeasibility for the above subproblem;
%% info.relgap = the relative duality gap for the above subproblem;
%% info.relkkt = the relative KKT residual for the above subproblem;
%% info.eta = the last value of eta;
%% info.nnzeros_x = the number of nonzero entries for x;
%% info.iterPPA = the total number of iterations for PPA;
%% info.iterSSN = the total number of iterations for SSN;
%% info.cntATmap = the total number of calls to ATmap;
%% info.cntAmap = the total number of calls to Amap;
%% info.time = total running time for dPPA_SSN;
%% info.ttime_cpu = total CPU time for dPPA_SSN;
%% runhist = a structure containing the history of the run;
%% dPPA_SSN:
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the Section 5 in the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%*************************************************************************
function [x,u,z,info,runhist] = dPPA_SSN(A,b,OPTIONS)

%% Input parameters
%%
maxiter_PPA = 200; % the maximum number of iterations for the PPA;
maxitersub_SSN = 200;  % the maximum number of iterations for the SSN;
etascale = 2; 
eta0 = 1; 
printSSN = 0;
printPPA = 0;
maxtime = 3600;

if isfield(OPTIONS,'m'), m = OPTIONS.m; else, m = length(b); end
if isfield(OPTIONS,'n'), n = OPTIONS.n; else, n = size(A,2); end
if isfield(OPTIONS,'kk1'), kk = OPTIONS.kk1; end
if isfield(OPTIONS,'lambda'), lambda = OPTIONS.lambda; end
if isfield(OPTIONS,'a'), a = OPTIONS.a; end
if isfield(OPTIONS,'rho'), rho = OPTIONS.rho; end
if isfield(OPTIONS,'sigma'), sigma = OPTIONS.sigma; end
if isfield(OPTIONS,'eta0'), eta0 = OPTIONS.eta0; end
if isfield(OPTIONS,'etascale'), etascale = OPTIONS.etascale; end
if isfield(OPTIONS,'obj_gap_mm'), obj_gap_mm = OPTIONS.obj_gap_mm;  end
if isfield(OPTIONS,'x0'), x0 = OPTIONS.x0; end
if isfield(OPTIONS,'u0'), u0 = OPTIONS.u0; end
if isfield(OPTIONS,'printPPA'), printPPA = OPTIONS.printPPA; end
if isfield(OPTIONS,'printSSN'), printSSN = OPTIONS.printSSN; end
if isfield(OPTIONS,'ATu'), ATu = OPTIONS.ATu; else, ATu = (u0'*A)'; end
if isfield(OPTIONS,'Ax'), Ax = OPTIONS.Ax;  end

par.m = m;
par.n = n;
par.kk = kk;
par.lambda = lambda;
par.rho = rho;
par.sigma = sigma;
par.x0 = x0;
par.a = a;
par.eta = eta0;
par.em = ones(m,1);

par.indexJ = logical(zeros(n,1));
par.AJ = [];
par.AJTAJ = [];
par.num_T = 0;

numSSN = 0;
breakPPA = 0;
tiny = 1e-10;
normb = norm(b);
cntATmap = 0;
cntAmap = 0;
epsilon = 1e-14;
%%
%% Amap and ATmap
%%
Amap = @(x) mexAx(A,x,0);
ATmap = @(y) (y'*A)';
%%
%% print the initial information
%%
tstart = clock;
tstart_cpu = cputime;
x_new = x0;
Ax0 = Ax;
Ax_snew = Ax0;

u_new = u0; z_new = Ax0-b;
Ax0_b = z_new;
dualobj_const = (par.rho/2)*norm(x0)^2+(par.sigma/2)*norm(Ax0_b)^2 + par.a'*Ax0_b;

w1 = x0 - (1/par.rho)*ATu;
xp = Proj_inf(w1,par.lambda/par.rho);
w2 = (1/par.sigma)*(u_new + a) + Ax0-b;
zp = Proj_dknorm(w2,par.kk/par.sigma,(1/par.sigma)*par.em);
gradphi = Amap(w1 - xp) - (w2 - zp) - b; % -gradphi

if obj_gap_mm > 5e-3
    update_rho = 2;
elseif obj_gap_mm > 1e-3
    update_rho = 3;
else
    update_rho = 4;
end

%%
%% =========================== Begin dPPA ==========================
%%
for iter = 1:maxiter_PPA
    par.iter = iter;
    w1_sub = w1; w2_sub = w2; %xp_sub = xp;  zp_sub = zp;
    u_snew = u_new; x_snew = x_new; z_snew = z_new; 
    gradPhi = gradphi; ATu_snew = ATu;
    Phi = -(par.rho/2)*norm(x_snew)^2 - (par.sigma/2)*norm(z_snew)^2 - u_snew'*b;
    parsub = par; normdelu = 0;
    
    %% ------------Begin SSN to compute: u_new -------------------
    for itersub = 1:maxitersub_SSN
        
        %% stopping critera: (C) and (D)
        normgradPhi = norm(gradPhi);
        subhist.gradPhi(itersub) = normgradPhi;
        const1 = 1/parsub.eta;
        theta = 1/(iter^(1.2)); 
        xi = min(1/(iter^(1.2)),0.5); 
        tolsub1 = const1*theta;
        if itersub == 1
            tolsub = min([0.1,tolsub1]);
        else
            tolsub2 = const1*xi*normdelu;
            tolsub = min([0.1,tolsub1,tolsub2]);
        end
        
        if ((normgradPhi <= tolsub) && (itersub > 1))  || (normgradPhi < 1e-9)
            if printSSN
                msg_SSN = 'good termination for SSN: ';
                fprintf('\n\t\t\t    %s  ',msg_SSN);
                fprintf('iterss=%2.0d,  gradPhi=%3.2e, tolsub=%3.2e, const1=%3.2e',...
                    itersub, normgradPhi, tolsub, const1);
            end
            break;
        end
        
        if (itersub > 100)
            ratio_gradphi = subhist.gradPhi(itersub-19:itersub)./subhist.gradPhi(itersub-20:itersub-1);
            if (min(ratio_gradphi) > 0.997) && (max(ratio_gradphi) < 1.003)
                if (printSSN);  fprintf('stagnate'); end
                break;
            end
        end
        
        %% Choose linsolver
        [flag_case,r_tmp,parsub] = Generatedash_Jacobi_MM(A,w1_sub,w2_sub,parsub);
        upnumber = 10000; r_indexJ = parsub.r_indexJ;
        subhist.r_indexJ(itersub) = r_indexJ;
        if (m <= upnumber)
            if (m <= 1000)
                linsolver = 'd_direct';
            elseif (r_indexJ <= max(0.01*n, upnumber))
                linsolver = 'd_direct';
            end
        end
        if (((r_indexJ <= m)||(r_indexJ < 4000)) && (r_indexJ <= upnumber))...
                || ((m >= 3000) && ((r_tmp <= 6500) || (r_indexJ <= 7000))) 
            linsolver = 'p_direct';
        end
              
        %% Compute Newton direction
        rhs = gradPhi;
        if strcmp(linsolver,'d_direct')
            if r_tmp == 0
                direc = parsub.eta*rhs;
            else
                W =  Cumpute_matrix_W_MM(flag_case,parsub);
                Hs = (1/parsub.rho)*(parsub.AJ*parsub.AJ') + (1/parsub.sigma)*W + (1/parsub.eta)*speye(m);
                if parsub.eta > 1e3
                    Hs = Hs + epsilon;
                end
                if m < 1500
                    direc = Hs\rhs;
                else
                    L = CholHess(Hs);
                    direc = linsysolvefun(L,rhs);
                end
            end
        elseif strcmp(linsolver, 'p_direct')
            if (parsub.r_indexJ == 0)
                [direc,~,~] = invDvec_MM(rhs,flag_case,parsub,0);
            else
                [invDrhs,T,parsub] = invDvec_MM(rhs,flag_case,parsub,1);
                if parsub.r_indexJ <= 1%1500
                    d_tmp = T\((invDrhs)'*parsub.AJ)';
                else
                    L = CholHess(T);
                    d_tmp = linsysolvefun(L,((invDrhs)'*parsub.AJ)');
                end
                [invDrhs_tmp,~,~] = invDvec_MM(parsub.AJ*d_tmp,flag_case,parsub,0);
                direc = invDrhs - invDrhs_tmp;
                parsub.num_T = parsub.num_T + 1;
            end
        end
        
        %% Strongly Wolfe search for stepsize: update u_snew
        steptol = 1e-7;
        [u_snew,x_snew,z_snew,alpha,iterstep,g0,w1_sub,w2_sub,Phi,maxiterfs,ATdu,normdelu]= findstep_new_MM...
            (parsub,ATmap,b,w1_sub,w2_sub,x_snew,z_snew,u_snew,u_new,direc,gradPhi,Phi,steptol);
        numSSN = numSSN + 1;
        cntATmap = cntATmap + 1;
        if alpha < tiny
            breakPPA =11; 
            msg_PPA = 'small steplength';
            break;
        end
        
        %% update -gradphi
        Ax_snew = Amap(x_snew); cntAmap = cntAmap + 1;
        gradPhi = Ax_snew - z_snew - b - (1/parsub.eta)*(u_snew - u_new);
        ATu_snew = ATu_snew + alpha*ATdu;        
        
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
         end
         
         runhist.phi(itersub) = -Phi;
         if (printSSN) && (rem(itersub,print_itersub)==0 || itersub == maxitersub_SSN)
            fprintf('\n\t\t\t\t [%2.0d]    [%3.2e   %3.2e]',itersub, normgradPhi, tolsub);
            fprintf('    flag=%1.1f     g0=%3.2e     lin=%2.0d  ', flag_case, -g0, solver);
            fprintf('   [%3.2e %2.0f] ',alpha,iterstep);
            %------------test----------------------------------------------
            fprintf('    [%3.0f  %3.0f  %3.0f  %3.0f] ', parsub.r_indexJ,...
                parsub.index_alpha,parsub.index_beta,parsub.index_gamma);
            if itersub > 1
                delphi = runhist.phi(itersub)-runhist.phi(itersub-1);
                fprintf('   delphi=%3.2e ', delphi);
            end
            %-----------end------------------------------------------------
            if maxiterfs == 1
                fprintf('$');
            end
         end
    end
    %% ------------------------- End SSN  -------------------------
    par.AJ = parsub.AJ; par.AJTAJ = parsub.AJTAJ; par.indexJ = parsub.indexJ;
    par.num_T = parsub.num_T;
    
    %% update u_new, x_new and z_new
    u_new = u_snew;
    x_new = x_snew;
    z_new = z_snew;
    Ax_new = Ax_snew;
    w1 = w1_sub; 
    w2 = w2_sub;
    
    %% stopping criteria
    Ax_b = Ax_new - b;
    def_Ax = Ax_new-Ax0; normsq_def_Ax = norm(def_Ax)^2;
    gradphi = Ax_b - z_new;
    normgradphi = norm(gradphi);
    
    gradphi_tmp = sort(abs(gradphi),'descend');
    err_ppa = sum(gradphi_tmp(1:par.kk)) + par.sigma*normgradphi^2;
    tol_ppa = (par.sigma/2)*normsq_def_Ax;
    
    % compute primfeas and relgap
    primfeas = normgradphi/(1 + normb);
    sort_Ax_b = sort(abs(Ax_b),'descend');
    primobj = sum(sort_Ax_b(1:kk)) + par.lambda*norm(x_new,1) - par.a'*def_Ax...
        + (par.rho/2)*norm(x_new-x0)^2 + (par.sigma/2)*normsq_def_Ax;
    dualobj = Phi + (1/(2*par.eta))*normdelu^2 + dualobj_const;
    relgap = abs(primobj-dualobj)/(1+abs(primobj)+abs(dualobj));
    
    % compute relkkt
    ATu = ATu_snew;
    x_tmp = par.rho*(x_new-x0) + ATu;
    z_tmp = par.a + par.sigma*(Ax0_b-z_new) + u_new;
    
    eta_u = primfeas;
    eta_x = norm(x_tmp + Proj_inf(x_new-x_tmp,par.lambda))/(1+norm(x_tmp)+norm(x_new));
    eta_z = norm(z_tmp  - Proj_dknorm(z_new + z_tmp,par.kk,par.em));
    relkkt = max([eta_u, eta_x, eta_z]);
    if (err_ppa <= tol_ppa) || (err_ppa < 5e-10)% || (relkkt < 1e-7)
        msg_PPA = 'dPPA converged';
        breakPPA = 1;
    end
    
    runhist.err_ppa(iter) = err_ppa;
    
    if (iter > 100)
        ratio_err_paa =  runhist.err_ppa(iter-9:iter)./ runhist.err_ppa(iter-10:iter-1);
        if (min(ratio_err_paa) > 0.997) && (max(ratio_err_paa) < 1.003)
            msg_PPA = 'stagnate PPA';
            breakPPA = 1111;
        end
    end
    
    %% Print results of dual PPA

    if printPPA == 1
        fprintf('\n\t\t {%3.0d}   {%3.2e %3.2e}    {%3.2e %3.2e}    %3.2e   {%3.2e %3.2e %3.2e}',...
            iter, err_ppa, tol_ppa, primfeas, relgap, relkkt, par.rho, par.sigma, par.eta);
    end
    
    %% termination
    ttime = etime(clock,tstart);
    if  (iter == maxiter_PPA)
        msg_PPA = ' maximum iteration reached';
        breakPPA = 100;
    end
    if  (ttime > maxtime)
        msg_PPA = ' maximum time reached';
        breakPPA = 1000;
    end
    if (breakPPA > 0) || (ttime > maxtime) || (iter == maxiter_PPA)
        if printPPA == 1
            fprintf('\n\t    breakPPA = %3.1f, %s, err_PPA = %3.2e',breakPPA,msg_PPA,err_ppa);
        end
        x = x_new; z = z_new; u = u_new;
        break;
    end
    
    %% Updata parameter
    if mod(iter,update_rho) == 0
        par.eta = min(etascale*par.eta,1e6);
    end
end
%%
%% =========================== End dPPA ============================
%%

%%
%% record history
%%
ttime = etime(clock,tstart);
ttime_cpu = cputime - tstart_cpu;

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

info.err_ppa = err_ppa;
info.tol_ppa = tol_ppa;
info.primfeas = primfeas;
info.relgap = relgap;
info.relkkt = relkkt;
info.eta = par.eta;
info.nnzeros_x = nnzeros_x;
info.iterPPA = iter;
info.iterSSN = numSSN;
info.cntATmap = cntATmap;
info.cntAmap = cntAmap + 1;

info.time = ttime;
info.ttime_cpu = ttime_cpu;
info.eta = par.eta;
runhist.x = x;
runhist.ATu = ATu;
runhist.Ax_b = Ax_b;
runhist.Ax = Ax_new;



