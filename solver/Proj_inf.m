function w = Proj_inf(w,r,b)
%% Compute the metric projection of the vector w on the ball
%% with center 0 and radius r in sense of ell_{inf} norm,
%% where the vector b is the weight of w
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


