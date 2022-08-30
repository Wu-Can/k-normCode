function psi_eps = psi_eps_fun(s,epsilon,flag)
%% If flag = 0, calculate the function value of psi;
%% if flag = 1, calculate the first-order derivative of psi;
%% if flag = 2, calculate the second-oeder derivative of psi;

if flag == 0
    psi_eps = (sqrt(4*epsilon^2 + s.^2) + s)/2; 
elseif flag == 1
    psi_eps = 0.5*(s./sqrt(4*epsilon^2 + s.^2) + 1);
elseif flag == 2
    psi_eps = 2*epsilon^2./((4*epsilon^2 + s.^2).^(1.5));
end