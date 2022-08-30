%% compute matrix*vector: (I_n+A^T*A)*x
% function Mx = Matvecnx(ATAmap,B,x)
% Mx = B^2*x + feval(ATAmap,x);
% end

function Mx = Matvecnx(ATAmap,x)
Mx = x + feval(ATAmap,x);
end