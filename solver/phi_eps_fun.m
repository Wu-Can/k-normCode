function phi_eps = phi_eps_fun(t,epsilon,flag)
%% If flag = 0, calculate the function value of phi;
%% if flag = 1, calculate the first-order derivative of phi;
%% if flag = 2, calculate the second-oeder derivative of phi;

if flag == 0
    phi_eps = sqrt(4*epsilon^2 + t.^2);
elseif flag == 1
    phi_eps = t./sqrt(4*epsilon^2 + t.^2);
elseif flag == 2
    phi_eps = 4*epsilon^2./((4*epsilon^2 + t.^2).^(1.5));
end