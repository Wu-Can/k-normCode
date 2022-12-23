%%=========================================================================
%% findstep:
%% find the step length alpha of SSN by the strong Wolfe line search 
%%
%% [u,x_snew,z_snew,alpha,iter,g0,w1_new,w2_new,Ly,maxiterfs,ATdu,normdelu]= ...
%%    findstep(par,ATmap,b,w1_sub,w2_sub,x_sub,z_sub,u0, u_old, du,gradphi,Ly0,tol,options)
%%
%% Input:
%% par.sigma = the parameter sigma in Algorithm 2 of the paper
%% par.tau = the parameter tau in Algorithm 2 of the paper
%% par.lambda = the parameter lambda in the model (9) of the paper
%% par.kk = the parameter k of k-norm
%% par.Im = the m-dimensional vector of all ones
%% ATmap = the mapping satisfying ATmap(y)= AT*y
%% b = the response vector
%% w1_sub = the vector omega_1 at u0 defined in Section 3.3 of the paper
%% w2_sub = the vector omega_2 at u0 defined in Section 3.3 of the paper
%% x_sub = w1_sub-Proj_inf(w1_sub,r0_1) with r0_1 = sigma*lambda
%% z_sub = w2_sub-Proj_dknorm(w2_sub,r0_2,sigma*Im) with r0_2 = k*sigma
%% u0 = the current iteration point u in Algorithm 2 of the paper
%% u_old = the initial iteration point u in Algorithm 2 of the paper
%% du = the Newton direction
%% gradphi = the right-hand term in equation (14) of the paper
%% Ly0 = the negative objective value of subproblem (11) at u0 
%%       in the paper
%% tol = a small positive number
%% Output:
%% u = the latest iteration point u in Algorithm 2 of the paper
%% x_snew = w1_new-Proj_inf(w1_new,r0_1)
%% z_snew = w2_new-Proj_dknorm(w2_new,r0_2,sigma*Im)
%% alpha = the step length
%% iter = the number of iterations of the strong Wolfe line search 
%% g0 = du'*gradphi
%% w1_new =the new vector omega_1 at u defined in Section 3.3 of the paper
%% w2_new = the new vector omega_2 at u defined in Section 3.3 of the paper
%% Ly = the negative objective value of subproblem (11) at u in the paper
%% maxiterfs = 1, if iter = maximum number of iterations
%%           = 0, otherwise
%% ATdu = AT*du
%% normdelu = norm(u-u_old)
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the Section 3 of the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%=========================================================================
function [u,x_snew,z_snew,alpha,iter,g0,w1_new,w2_new,Ly,maxiterfs,ATdu,normdelu]= ...
    findstep(par,ATmap,b,w1_sub,w2_sub,x_sub,z_sub,u0, u_old, du,gradphi,Ly0,tol,options)
if isfield(options,'stepop'); stepop = options.stepop; end
printlevel = 0;
maxit = ceil(log(1/(tol+eps))/log(2));
tol = 1e-6;
c1 = 1e-4; c2 = 0.9;
sigma = par.sigma; tau = par.tau; lambda = par.lambda;
kk = par.kk; Im = par.Im;
%%
g0 = du'*gradphi;
ATdu = feval(ATmap,du); 
tmp1 = sigma*ATdu; tmp2 = sigma*du;
bdu = b'*du; 
tmp3 = (tau/sigma)*norm(du)^2;
tmp4 = (tau/sigma)*(du'*(u0-u_old));
Ly = [];
normdelu = 0;
maxiterfs = 0;
r0_1 = sigma*lambda; r0_2 = kk*sigma;
if (g0 <= 0)
    alpha = 0; iter = 0;
    if (printlevel)
        fprintf('\n Need an ascent direction, %2.1e  ',g0);
    end
    u = u0; x_snew = x_sub; z_snew = z_sub;
    w1_new = w1_sub; w2_new = w2_sub;
    Ly = Ly0; normdelu =0;
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
    
    % update phi_new and gradphi_new
    w1_new = w1_sub-alpha*tmp1;
    xp_new = Proj_inf(w1_new,r0_1);
    w2_new = w2_sub+alpha*tmp2;
    zp_new = Proj_dknorm(w2_new,r0_2,sigma*Im); 
    x_snew = w1_new-xp_new;
    z_snew = w2_new-zp_new;
    
    galp = ATdu'*x_snew-du'*z_snew-bdu-alpha*tmp3-tmp4;
    
    if (iter==1)
        gLB = g0; gUB = galp;
        if (sign(gLB)*sign(gUB) > 0)
            if (printlevel); fprintf('|'); end
            normdelu = norm(u-u_old);
            Ly = -(1/(2*sigma))*(norm(x_snew)^2+norm(z_snew)^2)-b'*u-...
                (tau/(2*sigma))*normdelu^2;
            return;
        end
    end
    
    if (abs(galp) < c2*abs(g0))
        normdelu = norm(u-u_old);
         Ly = -(1/(2*sigma))*(norm(x_snew)^2+norm(z_snew)^2)-b'*u-...
                (tau/(2*sigma))*normdelu^2;
        if (Ly-Ly0-c1*alpha*g0 > -1e-8/max(1,abs(Ly0)))
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
if isempty(Ly)
    normdelu = norm(u-u_old);
    Ly = -(1/(2*sigma))*(norm(x_snew)^2+norm(z_snew)^2)-b'*u-...
                (tau/(2*sigma))*normdelu^2;
end
%%********************************************************************