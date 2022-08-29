%%*************************************************************************
%% psqmr:  preconditioned symmetric QMR with left (symmetric) preconditioner. 
%% Ax = b
%% b = rhs vector.
%% resnrm = norm of qmr-generated residual vector b-Ax. 
%%
%% This is for solving the linear system in ADMM_primal.m
%%*******************************************************************************
 function [x,Ax,resnrm,solve_ok] = psqmr_ADMM_new(Matvec,map,b,par,Ax,L) 
   
   N = length(b); 
   if ~exist('maxit') 
       maxit = max(100,sqrt(length(b))); % maximum iteration of PCG
   end
   if ~exist('L') 
       par.precond = 0;
       L = 1;
   end
   if isfield(par,'printlevel')
       printlevel = par.printlevel;
   else
       printlevel = 1;
   end
   if isfield(par,'tol')
       tol = par.tol;% the termination condition error in SSN
   else
       tol = 1e-6*norm(b); 
   end
   if isfield(par,'u0')
       x0 = par.u0;
   else
       x0 = zeros(N,1);
   end
   solve_ok = 1; 
   stagnate_check = 20; % ???
   miniter = 5;
   if isfield(par,'stagnate_check_psqmr')
      stagnate_check = par.stagnate_check_psqmr; 
   end
  if isfield(par,'minitpsqmr')
      miniter = par.minitpsqmr;
  end

%% compute r0
   x = x0; 
   if (norm(x) > 0) 
      Aq = Ax;
      r = b - Aq;
   else
      r = b;
   end
   err = norm(r); 
   resnrm(1) = err; 
   minres = err; 
%%
   q = precondfun(par,L,r); % z_0 = M^{-1}r
   tau_old  = norm(q);      
   rho_old  = r'*q; 
   theta_old = 0; 
   d = zeros(N,1); 
   res = r; 
   Ad = zeros(N,1);
%%      
%% main loop
%%
   tiny = 1e-30; 
   for iter = 1:maxit 
       Aq = feval(Matvec,map,q);
       sigma = q'*Aq; 
       if (abs(sigma) < tiny)
           solve_ok = 2;
           if (printlevel); fprintf('s1'); end
           Ax = feval(Matvec,map,x);
           break;
       else
          alpha = rho_old/sigma; % update stepsize alpha_k
          r = r - alpha*Aq;      % update r_k
       end
       u = precondfun(par,L,r);  % update z
       
       %% x = x + alpha*q, where alpha is stepsize, and q is a descend direction
       theta = norm(u)/tau_old; c = 1/sqrt(1+theta^2); 
       tau = tau_old*theta*c;
       gam = (c^2*theta_old^2); eta = (c^2*alpha); 
       d = gam*d + eta*q;
       x = x + d;                %% update x_k
       %%----- stopping conditions ----
       Ad = gam*Ad + eta*Aq;
       res = res - Ad; %update r
       err = norm(res); 
       resnrm(iter+1) = err; 
       if (err < minres); minres = err; end
       if (err < tol) && (iter > miniter)
           Ax = feval(Matvec,map,x);
           break; 
       end  
       % check the stagnation point
       if (iter > stagnate_check) && (iter > 10)
           ratio = resnrm(iter-9:iter+1)./resnrm(iter-10:iter); 
           if (min(ratio) > 0.997) && (max(ratio) < 1.003)
               solve_ok = -1; 
               if (printlevel);  fprintf('s'); end
               Ax = feval(Matvec,map,x);
               break;
           end
       end
       %%----------------------------- 
       if (abs(rho_old) < tiny)
          solve_ok = 2;
          if (printlevel); fprintf('s2'); end
          Ax = feval(Matvec,map,x);
          break;
       else
          rho  = r'*u; 
          beta = rho/rho_old; % update beta
          q = u + beta*q; % update direction
       end
       rho_old = rho; 
       tau_old = tau; 
       theta_old = theta; 
   end
   if (iter == maxit)
       solve_ok = -2; 
       Ax = feval(Matvec,map,x);
   end
   if (solve_ok ~= -1)
       if (printlevel)
           fprintf(' '); 
       end
   end
   
%%************************************************************************
%%************************************************************************
   function  q = precondfun(par,L,r)
   precond = 0; 
   if isfield(par,'precond')
      precond = par.precond; 
   end
   
   if (precond == 0)
      q = r;
   elseif (precond == 1)
      q = L.invdiagM.*r;
   end
%%************************************************************************
