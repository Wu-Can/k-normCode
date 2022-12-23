%%=========================================================================
%% psqmr_knorm_N_ALM:
%% preconditioned symmetric QMR with left (symmetric) preconditioner
%% for solving the Newton linear system in NALM.m
%% (E)                Ax = b
%%
%% [x,resnrm,solve_ok] = psqmr_knorm_N_ALM(matvecfname,C,b,par,L,tol,
%%                       maxit,printlevel) 
%%
%% Input: 
%% matvecfname = 'matvec_N_ALM'
%% C = part of matrix A (C*C^T = AJ*AJ^T + W)
%% b = the vector in (E)
%% par.precond = 0, does not use preconditioner; otherwise, use it
%% par.d = initial point x0 
%% par.stagnate_check_psqmr = the minimum number of iterations to 
%%                             check the stagnation point
%% par.minitpsqmr = the minimum number of iterations for the psqmr
%% L = the information of the preconditioner
%% tol = tolerance for the stopping criterion
%% maxit = the maximum number of iterations
%% printlevel = 0, does not print the status of the solution of (E);
%%                otherwise, print it
%% Output:
%% x = the output solution of (E)
%% resnrm = the vector of residuals norm(b-Ax) for each iteration of 
%%          the psqmr 
%% solve_ok = the information of finial status
%%=========================================================================
 function [x,resnrm,solve_ok] = psqmr_knorm_N_ALM(matvecfname,C,b,...
       par,L,tol,maxit,printlevel) 
   
   N = length(b); 
   if ~exist('maxit') 
       maxit = max(50,sqrt(length(b))); % maximum iteration of PCG
   end
   if ~exist('printlevel')
       printlevel = 0; 
   end
   if ~exist('tol')
       tol = 1e-6*norm(b); % the termination condition error in SSN
   end
   if ~exist('L') 
       par.precond = 0; 
   end
  
   x0 = zeros(N,1);
   if isfield(par,'d')
       x0 = par.d;
   end
   solve_ok = 1; 
   stagnate_check = 20; 
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
      Aq = feval(matvecfname,x,C,par);
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
        Aq = feval(matvecfname,q,C,par);
       sigma = q'*Aq; 
       if (abs(sigma) < tiny)
          solve_ok = 2;
          if (printlevel); fprintf('s1'); end
         break;
       else
          alpha = rho_old/sigma; % update step size alpha_k
          r = r - alpha*Aq;      
       end
       u = precondfun(par,L,r);  
       
       %% x = x + d
       theta = norm(u)/tau_old; c = 1/sqrt(1+theta^2); 
       tau = tau_old*theta*c;
       gam = (c^2*theta_old^2); eta = (c^2*alpha); 
       d = gam*d + eta*q;
       x = x + d;               
       
       %% stopping conditions
       Ad = gam*Ad + eta*Aq;
       res = res - Ad; 
       err = norm(res); 
       resnrm(iter+1) = err; 
       if (err < minres); minres = err; end
       if (err < tol) && (iter > miniter)
           if (printlevel);  fprintf('s3'); end
           break; 
       end  
       %% check the stagnation point
       if (iter > stagnate_check) && (iter > 10) 
           ratio = resnrm(iter-9:iter+1)./resnrm(iter-10:iter); 
           if (min(ratio) > 0.997) && (max(ratio) < 1.003)
               solve_ok = -1; 
               if (printlevel);  fprintf('s'); end
               break;
           end
       end
       %%----------------------------- 
       if (abs(rho_old) < tiny)
          solve_ok = 2;
          if (printlevel); fprintf('s2'); end
          break;
       else
          rho  = r'*u; 
          beta = rho/rho_old;
          q = u + beta*q;
       end
       rho_old = rho; 
       tau_old = tau; 
       theta_old = theta; 
   end
   if (iter == maxit)
       solve_ok = -2; 
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
   if isempty(L.invdiagM)
      precond = 0; 
   end  
   if (precond == 0) 
      q = r;
   elseif (precond == 1)
      q = L.invdiagM.*r;
   end
%%************************************************************************
