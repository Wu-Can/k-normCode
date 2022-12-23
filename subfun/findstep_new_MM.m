%%=========================================================================
%% findstep_new_MM:
%% find the step length alpha of SSN by the strong Wolfe line search
%% 
%% [u,x_snew,z_snew,alpha,iter,g0,w1_new,w2_new,Lu,maxiterfs,ATdu,normdelu]= ...
%%    findstep_new_MM(par,ATmap,b,w1_sub,w2_sub,x_sub,z_sub,u0,u_old, du,gradPhi,Lu0,tol)
%%
%% Input:
%% par.rho = the parameter rho in (24) of the paper
%% par.sigma = the parameter sigma in (24) of the paper
%% par.eta = the parameter eta of the PPA
%% par.lambda = the parameter lambda in (23) of the paper
%% par.kk = the parameter k_1 in f_1 of (23) in the paper
%% par.em = the m-dimensional vector of all ones
%% ATmap = the mapping satisfying ATmap(y)= AT*y
%% b = the response vector
%% w1_sub = a n-dimensional vector
%% w2_sub = a m-dimensional vector
%% x_sub = w1_sub-Proj_inf(w1_sub,r0_1) with r0_1 = lambda/rho
%% z_sub = w2_sub-Proj_dknorm(w2_sub,r0_2,(1/sigma)*em) 
%%         with r0_2 = k_1/sigma
%% u0 = the current iteration point u in the SSN
%% u_old = the initial iteration point u in the SSN
%% du = the Newton direction
%% gradPhi =  the right-hand term of the Newton system
%% Lu0 = the negative objective value of the subproblem of the PPA
%%       at u0
%% tol = a small positive number
%% Output:
%% u = the latest iteration point u in the SSN
%% x_snew = w1_new-Proj_inf(w1_new,r0_1) 
%% z_snew = w2_new-Proj_dknorm(w2_new,r0_2,(1/sigma)*em)
%% alpha = the step length 
%% iter = the number of iteration of the strong Wolfe line search 
%% g0 = du'*gradPhi
%% w1_new = a new n-dimensional vector w1
%% w2_new = a new m-dimensional vector w2
%% Lu = the latest negative objective value of the subproblem of the
%%      PPA at u
%% maxiterfs = 1, if iter = maximum number of iterations
%%           = 0, otherwise
%% ATdu = AT*du
%% normdelu = norm(u-u_old)
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the Section 5 of the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%=========================================================================
function [u,x_snew,z_snew,alpha,iter,g0,w1_new,w2_new,Lu,maxiterfs,ATdu,normdelu]= ...
    findstep_new_MM(par,ATmap,b,w1_sub,w2_sub,x_sub,z_sub,u0, u_old, du,gradPhi,Lu0,tol)
printlevel = 0;
maxit = ceil(log(1/(tol+eps))/log(2));
c1 = 1e-4; c2 = 0.9;
rho = par.rho; sigma = par.sigma; eta = par.eta;
lambda = par.lambda;
kk = par.kk;
em = par.em;
%%
g0 = du'*gradPhi;
ATdu = feval(ATmap,du);
tmp1 = (1/rho)*ATdu; tmp2 = (1/sigma)*du;
bdu = b'*du;
tmp3 = (1/eta)*norm(du)^2;
tmp4 = (1/eta)*(du'*(u0-u_old));
Lu = [];
normdelu = 0;
maxiterfs = 0;
r0_1 = lambda/rho; r0_2 = kk/sigma;
if (g0 <= 0)
    alpha = 0; iter = 0;
    if (printlevel)
        fprintf('\n Need an ascent direction, %2.1e  ',g0);
    end
    u = u0; x_snew = x_sub; z_snew = z_sub;
    w1_new = w1_sub; w2_new = w2_sub;
    Lu = Lu0; normdelu =0;
    return;
end
%%
alpha = 1; alpconst = 0.5;
for iter = 1:maxit
    if (iter==1)
        alpha = 1; LB = 0; UB = 1;
    else
        alpha = alpconst*(LB+UB);
    end
    u = u0 + alpha*du;
    
    % update Phi_new and gradPhi_new
    w1_new = w1_sub-alpha*tmp1;
    xp_new = Proj_inf(w1_new,r0_1);
    w2_new = w2_sub+alpha*tmp2;
    zp_new = Proj_dknorm(w2_new,r0_2,(1/sigma)*em);
    x_snew = w1_new-xp_new;
    z_snew = w2_new-zp_new;
    
    galp = ATdu'*x_snew-du'*z_snew-bdu-alpha*tmp3-tmp4;
    
    if (iter==1)
        gLB = g0; gUB = galp;
        if (sign(gLB)*sign(gUB) > 0)
            if (printlevel); fprintf('|'); end
            normdelu = norm(u-u_old);
            Lu = -(rho/2)*norm(x_snew)^2 - (sigma/2)*norm(z_snew)^2 - b'*u-...
                (1/(2*eta))*normdelu^2;
            return;
        end
    end
    
    if (abs(galp) < c2*abs(g0))
        normdelu = norm(u-u_old);
        Lu = -(rho/2)*norm(x_snew)^2 - (sigma/2)*norm(z_snew)^2 - b'*u-...
            (1/(2*eta))*normdelu^2;
        if (Lu-Lu0-c1*alpha*g0 > -1e-8/max(1,abs(Lu0)))
            % && ((stepop==1) || (stepop==2 && abs(galp)<tol))
            if (printlevel); fprintf(':'); end
            return;
        end
    end
    
    if (sign(galp)*sign(gUB) < 0)
        LB = alpha; gLB = galp;
    elseif (sign(galp)*sign(gLB) < 0)
        UB = alpha; gUB = galp;
    end
end
if (iter == maxit)
    maxiterfs = 1;
end
if (printlevel); fprintf('m'); end
if isempty(Lu)
    normdelu = norm(u-u_old);
    Lu = -(rho/2)*norm(x_snew)^2 - (sigma/2)*norm(z_snew)^2 - b'*u-...
        (1/(2*eta))*normdelu^2;
end
%%********************************************************************