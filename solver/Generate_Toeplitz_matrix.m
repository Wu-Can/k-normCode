function T = Generate_Toeplitz_matrix(a,n)
%% generate n*n Toeplitz covariance matrix T with T(i,j)=a^(|i-j|)
a_mul = 1;
for j = 1:n
    a_mul = a_mul*a;
    if a_mul < eps
        nnz_x = j;
        break;
    end
    nnz_x = j;
end
nnz_x = min(nnz_x,n-1);

x_tmp = zeros(nnz_x,1);
for i =1:nnz_x
    x_tmp(i) = a^i;
end

x = zeros(n,1);
if n <= nnz_x
    x = x_tmp(1:n);
else
    x(1:nnz_x) = x_tmp;
end

T_tmp = sparse(n,n);

for jj = 1:n-1
    t = min(n,nnz_x+jj);
    T_tmp(jj+1:t,jj) = x(1:t-jj);
end
T = T_tmp + speye(n) + T_tmp';

% T_2 = a.^toeplitz(0:n-1);
% error = norm(T-T_2);
end

