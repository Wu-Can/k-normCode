%%*******************************************************************************
%% qp_Gurobi:
%% Run the barrier method in Gurobi for solving the nonconvex 
%% truncated CVaR-based sparse linear regression under QP form:
%% minimize    k*c + sum(u) + lambda*sum(v) - < AT*a, x - x_old> 
%% {x,v,u,c}   + (rho/2)||x - x_old||^2 + (sigma/2)||A*x - A*x_old||^2
%% subject to             A*x - u - c*e <= b
%%                        A*x + u + c*e >= b
%%                                x + v >= 0
%%                                x - v <= 0
%%                                    u >= 0
%% where x, v in R^n, u in R^m, c in R are varibles and e:=ones(m,1)
%%*******************************************************************************
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



