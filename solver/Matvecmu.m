%% compute matrix*vector: (I_m+A*A')*x
function Mu = Matvecmu(AATmap,u)
Mu = u + feval(AATmap,u);
end

