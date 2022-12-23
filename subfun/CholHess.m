%%=========================================================================
%% CholHess:
%% Cholesky factorization of a symmetric positive definite matrix 
%% Hess such that Hess = L.R'*L.R with an upper triangular L.R         
%%
%% L = CholHess(Hess)
%%
%% Input:
%% Hess = a symmetric positive definite matrix
%% Output:
%% L.R = the upper triangular satisfying Hess = L.R'*L.R
%% L.Rt = the transpose of L.R 
%% L.flag = 0, Hess is a symmetric positive definite matrix and 
%%             the factorization was successful; otherwise not.
%% L.perm = the vector satisfying Hess(L.perm,L.perm) = L.R'*L.R 
%%          if Hess is sparse; otherwise, L.perm = [1:m], where m 
%%          is the dimension of the symmetric matrix Hess   
%% L.matfct_options = 'spCholsky_matlab', if Hess is sparse
%%                  = 'Cholesky', if Hess is not sparse
%%=========================================================================
function L = CholHess(Hess)
m = size(Hess,2);

if (nnz(Hess)<0.2*m*m)
    use_spchol = 1;
else
    use_spchol = 0;
end

if (use_spchol)
    [L.R, L.flag, L.perm] = chol(sparse(Hess),'vector'); %;
    L.Rt = L.R';
    L.matfct_options = 'spCholsky_matlab';
else
    if issparse(Hess)
        Hess = full(Hess);
    end
    L.matfct_options = 'Cholesky';
    L.perm = [1:m];
    [L.R, L.flag] = chol(Hess); % 
    L.Rt = L.R';
end

end
