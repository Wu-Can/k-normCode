%%=========================================================================
%% Matvecmu:
%% compute matrix*vector: (I_m+A*A')*u
%%
%% Mu = Matvecmu(AATmap,u)
%%
%% Input:
%% AATmap = a matrix-vector product mapping
%% u = a vector
%% Output:
%% Mu = (I_m+A*A')*u
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%=========================================================================
function Mu = Matvecmu(AATmap,u)
Mu = u + feval(AATmap,u);
end

