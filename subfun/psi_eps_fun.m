%%=========================================================================
%% psi_eps_fun:
%% calculate the function value, the first-order derivative or
%% the second-order derivative of the function
%%           psi_eps(t) = (sqrt(4*epsilon^2 + s.^2) + s)/2
%%
%% psi_eps = psi_eps_fun(s,epsilon,flag)
%%
%% Input: 
%% s = the variable of psi_eps
%% epsilon = the parameter in psi_eps
%% flag = 0, calculate the function value of psi_eps
%%      = 1, calculate the first-order derivative of psi_eps
%%      = 2, calculate the second-order derivative of psi_eps
%% Output:
%% psi_eps = the function value of psi_eps, if flag = 0
%%         = the first-order derivative of psi_eps, if flag = 1
%%         = the second-order derivative of psi_eps, if flag = 2
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the Appendix E in supplementary 
%% materials of the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%=========================================================================
function psi_eps = psi_eps_fun(s,epsilon,flag)
if flag == 0
    psi_eps = (sqrt(4*epsilon^2 + s.^2) + s)/2; 
elseif flag == 1
    psi_eps = 0.5*(s./sqrt(4*epsilon^2 + s.^2) + 1);
elseif flag == 2
    psi_eps = 2*epsilon^2./((4*epsilon^2 + s.^2).^(1.5));
end
