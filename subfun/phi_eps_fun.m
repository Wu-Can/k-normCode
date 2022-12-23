%%=========================================================================
%% phi_eps_fun:
%% calculate the function value, the first-order derivative or
%% the second-order derivative of the function
%%          phi_eps(t):=sqrt(4*epsilon^2+t.^2)
%%
%% phi_eps = phi_eps_fun(t,epsilon,flag)
%%
%% Input: 
%% t = the variable of phi_eps
%% epsilon = the parameter in phi_eps
%% flag = 0, calculate the function value of phi_eps
%%      = 1, calculate the first-order derivative of phi_eps
%%      = 2, calculate the second-order derivative of phi_eps
%% Output:
%% phi_eps = the function value of phi_eps, if flag = 0
%%         = the first-order derivative of phi_eps, if flag = 1
%%         = the second-order derivative of phi_eps, if flag = 2
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the Appendix E in supplementary 
%% materials of the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%=========================================================================
function phi_eps = phi_eps_fun(t,epsilon,flag)
if flag == 0
    phi_eps = sqrt(4*epsilon^2 + t.^2);
elseif flag == 1
    phi_eps = t./sqrt(4*epsilon^2 + t.^2);
elseif flag == 2
    phi_eps = 4*epsilon^2./((4*epsilon^2 + t.^2).^(1.5));
end
