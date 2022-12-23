%%*************************************************************************
%% k_norm_lp:
%% Run the barrier method in Gurobi for solving the convex 
%% CVaR-based sparse linear regression under LP form:
%%  (P)   minimize_{x,v,u,c}     k*c + sum(u) + lambda*sum(v)
%%        subject to             A*x - u - c*e <= b
%%                               A*x + u + c*e >= b
%%                               x + v >= 0
%%                               x - v <= 0
%%                               u >= 0
%% where x, v in R^n, u in R^m, c in R are varibles and e:=ones(m,1)
%%
%% [sol,fval,exitflag,output] = Gurobi_knorm_lp(A,b,OPTIONS)
%%
%% Input: 
%% A, b = matrix A and vector b in (P)
%% OPTIONS.m = the sample size
%% OPTIONS.n = the feature size
%% OPTIONS.kk = the value of parameter k in (P)
%% OPTIONS.lambda = the value of parameter lambda in (P)
%% Output:
%% sol.x = the output solution x in (P)
%% fval = the output objective value of (P)
%% exitflag = -3, problem is unbounded
%%          = -2, no feasible point found
%%          = 0, maximum number of iterations reached
%%          = 1, converged to a solution
%% output.baritercount = the number of iterations of the barrier 
%%                       method
%% output.time = total running time
%% output.constrviolation = maximum violation for constraints and
%%                          bounds
%%*************************************************************************
function [sol,fval,exitflag,output] = Gurobi_knorm_lp(A,b,OPTIONS)
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
[sol,fval,exitflag,output] = solve(prob, 'Options', options);

end