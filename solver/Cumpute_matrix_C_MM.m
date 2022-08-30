function C =  Cumpute_matrix_C_MM(flag_case,options)
%% This function is to calculate the m*m matrix C such that C*C' + (1/eta)*Im = Hessian,
%% where Hessian := (1/rho)AJ*AJT + (1/sigma)W + (1/eta)Im.

AJ = options.AJ;
index_alpha = options.index_alpha;
index_beta = options.index_beta;
Palpha = options.Palpha;
Pbeta = options.Pbeta;
Pgamma = options.Pgamma;
m = options.m;
rho = options.rho;
sigma = options.sigma;

num_1 = 1/sqrt(rho);
num_2 = 1/sqrt(sigma);

switch(flag_case)
    case 1
    C = num_1*AJ;
    case {2.1, 2.2}
    C = [num_1*AJ, num_2*(Palpha)'];
    case {3.1, 3.2}
    C = [num_1*AJ, num_2*(1/sqrt(index_alpha+index_beta))*(ones(1,index_alpha+index_beta)*[Palpha;Pbeta])', num_2*(Pgamma)'];
    case {4.1, 4.2}
    C = [num_1*AJ, num_2*(Palpha)', num_2*(1/sqrt(index_beta))*(ones(1,index_beta)*Pbeta)', num_2*(Pgamma)'];
    case 5
    C = [num_1*AJ,num_2*speye(m)];
end

end


