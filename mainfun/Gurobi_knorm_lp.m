%%*****************************************************************************
%% k_norm_lp:
%% Run the barrier method in Gurobi for solving the convex 
%% CVaR-based sparse linear regression under LP form:
%%     minimize_{x,v,u,c}     k*c + sum(u) + lambda*sum(v)
%%     subject to             A*x - u - c*e <= b
%%                            A*x + u + c*e >= b
%%                            x + v >= 0
%%                            x - v <= 0
%%                            u >= 0
%% where x, v in R^n, u in R^m, c in R are varibles and e:=ones(m,1).
%%*****************************************************************************
function [sol,fval,exitflag,output,lambda] = Gurobi_knorm_lp(A,b,OPTIONS)
m = OPTIONS.m;
n = OPTIONS.n;
k = OPTIONS.kk;
lambda = OPTIONS.lambda;

x = optimvar('x',n);
v = optimvar('v',n);
u = optimvar('u',m,1,'LowerBound',0);
c = optimvar('c');
c_tmp = c.*ones(m,1);

prob = optimproblem;
prob.Objective = k*c + sum(u)+lambda*sum(v);
prob.Constraints.cons1 = A*x-u-c_tmp <= b;
prob.Constraints.cons2 = A*x+u+c_tmp >= b;
prob.Constraints.cons3 = x+v >= 0;
prob.Constraints.cons4 = x-v <= 0;

options = optimoptions('linprog');
[sol,fval,exitflag,output,lambda] = solve(prob, 'Options', options);

end