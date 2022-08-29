function L = CholHess(Hess)
%% Cholesky Fractorization of Hess, i.e., Hess = L*L'
m = size(Hess,2);

if (nnz(Hess)<0.2*m*m)
    use_spchol = 1;
else
    use_spchol = 0;
end

if (use_spchol)
    [L.R, L.flag, L.perm] = chol(sparse(Hess),'vector'); %Hess(L.perm,L.perm)=R'*R;
    L.Rt = L.R';
    L.matfct_options = 'spCholsky_matlab';
else
    if issparse(Hess)
        Hess = full(Hess);
    end
    L.matfct_options = 'Cholesky';
    L.perm = [1:m];
    [L.R, L.flag] = chol(Hess); % A = L.R'*L.R; L.flag=0 when A is positive definite,otherwise L.flag>0
    L.Rt = L.R';
end

end