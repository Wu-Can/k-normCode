%%*************************************************************************
%% linsysolvefun: Solve H*x = b
%%
%% x = linsysolvefun(L,b)
%% where L contains the triangular factors of H.
%% 
%% SDPNAL: 
%% Copyright (c) 2008 by
%% Xinyuan Zhao, Defeng Sun, and Kim-Chuan Toh 
%%*************************************************************************
 
function x = linsysolvefun(L,b)
 x = zeros(size(b)); 
 if strcmp(L.matfct_options,'Cholesky')
     x(L.perm) = L.R \ (b(L.perm)' / L.R)'; % H=R'*R, x=R^{-1}(R^{-1})'b=R^{-1}(b'*R^{-1})'
 elseif strcmp(L.matfct_options,'spCholsky_matlab')
     x(L.perm) = mexbwsolve(L.Rt,mexfwsolve(L.R,b(L.perm))); 
 elseif strcmp(L.matfct_options,'lu')
     x(:) = L.u \ (L.l \ b(L.perm));
 elseif strcmp(L.matfct_options,'splu')     
     x(L.perm) = L.q*( L.u \ (L.l \ (L.p*b(L.perm))));
 end
end
%%*************************************************************************
