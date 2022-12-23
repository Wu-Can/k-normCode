%%=========================================================================
%% Matvecnx:
%% compute matrix*vector: (I_n+A^T*A)*x
%%
%% Mx = Matvecnx(ATAmap,x)
%%
%% Input:
%% ATAmap = a matrix-vector product mapping
%% x = a vector
%% Output:
%% Mx = (I_n+A^T*A)*x
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%=========================================================================
function Mx = Matvecnx(ATAmap,x)
Mx = x + feval(ATAmap,x);
end