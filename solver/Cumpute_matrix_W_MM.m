function W =  Cumpute_matrix_W_MM(flag_case,options)
%% This function is to calculate an elemrnt in
%% the Clarke generalized Jacobian of Prox_{(1/sigma)*||\cdot||_{k_1}}(\cdot)

index_alpha = options.index_alpha;
index_beta = options.index_beta;
index_gamma = options.index_gamma;
Palpha = options.Palpha;
Pbeta = options.Pbeta;
Pgamma = options.Pgamma;
m = options.m;

switch(flag_case)
    case 1
        W = sparse(m,m);
    case {2.1, 2.2}
        W = (Palpha)'*Palpha;
    case {3.1, 3.2}
        if index_gamma == 0
            eTP = ones(1,index_alpha+index_beta)*[Palpha;Pbeta];
            W = 1/(index_alpha+index_beta)*(eTP'*eTP);
        else
            eTP = ones(1,index_alpha+index_beta)*[Palpha;Pbeta];
            W = 1/(index_alpha+index_beta)*(eTP'*eTP) + Pgamma'*Pgamma;
        end
    case {4.1, 4.2}
        if index_gamma == 0
            eTP = ones(1,index_beta)*Pbeta;
            W = Palpha'*Palpha + (1/index_beta)*(eTP'*eTP);
        else
            eTP = ones(1,index_beta)*Pbeta;
            W = Palpha'*Palpha + (1/index_beta)*(eTP'*eTP) + Pgamma'*Pgamma;
        end
    case 5
        W = speye(m);
end

end


