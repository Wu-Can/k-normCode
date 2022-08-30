function Hessianv = matvec_N_ALM(d,C,options)
%%  Hess*d =-rhs,
%% This function is to generate the mapping Hess*d
sigma = options.sigma;
tau = options.tau;
    Cmap = @(x) mexAx(C,x,0);
    CTmap = @(x) mexAx(C,x,1);
    Hessianv = sigma*Cmap(CTmap(d))+(tau/sigma)*d; 
end