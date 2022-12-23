%%*************************************************************************
%% qp_Gurobi:
%% Run the barrier method in Gurobi for solving the nonconvex 
%% truncated CVaR-based sparse linear regression under QP form:
%% minimize    k*c + sum(u) + lambda*sum(v) - < AT*a, x - x_old> 
%% {x,v,u,c}   + (rho/2)||x - x_old||^2 + (sigma/2)||A*x - A*x_old||^2
%% subject to             A*x - u - c*e <= b
%% (P)                    A*x + u + c*e >= b
%%                                x + v >= 0
%%                                x - v <= 0
%%                                    u >= 0
%% where x, v in R^n, u in R^m, c in R are varibles and e:=ones(m,1)
%%
%% [x, fval, exitflag, output] = qp_Gurobi(A,b,options)
%%
%% Input:
%% A, b = the m*n-dimensional design matrix A and the m-dimensional 
%%        response vector b in (P)
%% options.m = sample size
%% options.n = feature size
%% options.kk1 = the value of parameter k in (P)
%% options.lambda = the value of parameter lambda in (P)
%% options.sigma = the value of parameter sigma in (P)
%% options.a = the m-dimensional vector a in (P)
%% options.ATA = the n*n-dimensional matrix AT*A
%% options.x0 = the n-dimensional vector x_old in (P)
%% Output:
%% x = the output solution of (P)
%% fval = the output objective value of (P)
%% exitflag = -3, problem is unbounded
%%          = -2, no feasible point found
%%          = 0, maximum number of iterations reached
%%          = 1, converged to a solution
%% output.message = message from Gurobi
%% output.time = total running time 
%% output.iterBar = the number of iterations of the barrier method
%% output.nnzeros_x = the number of nonzero elements of x
%%*************************************************************************
function [x, fval, exitflag, output] = qp_Gurobi(A,b,options)

m = options.m;
n = options.n;
k = options.kk1;
lambda = options.lambda;
rho = options.rho;
sigma = options.sigma;

a = options.a;
ATA = options.ATA;
x_old = options.x0;

ATAxv = ATA*x_old;
ATa = A'*a;

I_ATA = (rho/2)*speye(n) + (sigma/2)*ATA;

model.Q = sparse(2*n+m+1,2*n+m+1);
model.Q(1:n,1:n) = I_ATA;
model.A = sparse([A, sparse(m,n), speye(m), ones(m,1);
                 -A, sparse(m,n), speye(m), ones(m,1);
          -speye(n), speye(n), sparse(n,m), sparse(n,1);
           speye(n), speye(n), sparse(n,m), sparse(n,1)]);
model.obj = [-(rho*x_old' + sigma*ATAxv' + ATa'), lambda*ones(1,n), ones(1,m),k];
model.rhs = [b', -b', ones(1,n), ones(1,n)];
model.lb = [-inf*ones(1,n), -inf*ones(1,n), zeros(1,m), -inf];
model.sense = '>';

%% parameter setting
params.TimeLimit = 7200; % Time limit; default: infinity
params.BarConvTol = 1e-6; % Barrier convergence tolerance; default: 1e-8 (relative gap)
params.FeasibilityTol = 1e-6; % primal feasibility tol
params.OptimalityTol = 1e-6;  % dual feasiblility tol
params.TuneOutput = 3; % Tuning output level (2+detailed solver output);
params.Method = 2; % barrier method; % default: -1 (automatic)
params.Crossover = 0;

gurobi_write(model, 'qp.lp');
results = gurobi(model,params);

if strcmp(results.status,'OPTIMAL')
    exitflag = 1; % converged to a solution
elseif strcmp(results.status,'UNBOUNDED')
    exitflag = -3; % problem is unbounded
elseif strcmp(results.status,'ITERATION_LIMIT')
   exitflag = 0; % maximum number of iterations reached
else
    exitflag = -2; % no feasible point found
end

if ~(exitflag == 1)
    output.message = results.status;
    fprintf('Message from Gurobi: %s',output.message);
end

output.message = results.status;
output.time = results.runtime;
output.iterBar = results.baritercount;
x = [];
if isfield(results,'x')
    x = results.x(1:n);
    %----------------------the number of nonzeros in x--------------------
    sortx = sort(abs(x),'descend');
    normx1 = 0.999*norm(x,1);
    tmpidex = find(cumsum(sortx) > normx1);
    if isempty(tmpidex)
        nnzeros_x = 0;
    else
        nnzeros_x = tmpidex(1);
    end
    output.nnzeros_x = nnzeros_x;
end
fval = [];
if isfield(results,'objval')
    fval = results.objval;
end



