%%=========================================================================
%% Proj_inf:
%% compute the (weighted) projection of the vector w on the ball
%% with center 0 and radius r in the sense of ell_{inf} norm
%% 
%% w = Proj_inf(w,r,b)
%%
%% Input: 
%% w = a vector 
%% r = a non-negative constant as the radius
%% b = the weight of w
%% Output:
%% w = the (weighted) projection of the input vector w on the ball 
%%     with center 0 and radius r in the sense of ell_{inf} norm
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%=========================================================================
function w = Proj_inf(w,r,b)
k = length(r);
if k ~= 1
    fprintf('\n Error: the dimension of radius %2d should be 1',k)
    return;
end
if nargin == 2
    w = max(-r,min(r,w)); % better than  w(w > r) = r; w(w < -r) = -r;   
else
    w = max(-r*b,min(r*b,w));
end


