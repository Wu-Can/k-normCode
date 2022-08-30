function Hessianv = matvec_org(v,C,flag_case,flag_small,options)
%%  Hess*d =-rhs, 
%% generate mapping Hess*d
%% parameter : sig,tau,lam, k is the number in k-norm, and r is the radius, tol is a residual.

beta = options.beta;
sig = options.sig;
tau = options.tau;

Cmap = @(x) C*x;
CCTmap = @(x) Cmap((x'*C)');
CTCmap = @(x) ((Cmap(x))'*C)';

if flag_case == 5
    if flag_small
       Hessianv = sig*CCTmap(v) + (sig/tau)*v;    
   else
       Hessianv = CTCmap(v) + (1/tau)*v;       
    end 
elseif flag_case > 0
     if flag_small
       Hessianv = sig*CCTmap(v)+ (1/beta)*v;     
   else
       Hessianv = CTCmap(v)+(1/(sig*beta))*v;       
     end
else
     fprintf ('There not exists generalized Hessian, since flagcase = %3.0f', flag_case);
    return;
end    

end
    
    

