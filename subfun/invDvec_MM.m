%%===========================================================================
%% invDvec_MM:
%% calculate the vector D^{-1}*rhs 
%% and/or the matrix T = rho*I_s + AJ^{T}*D^{-1}*AJ
%%
%% [invDrhs,T,options] = invDvec_MM(rhs,flag_case,options,computeT)
%%
%% Input:
%% rhs = a vector 
%% flag_case = a positive number depending on the form of the matrix D
%% options.sigma = the value of parameter sigma 
%% options.eta = the value of parameter eta
%% options.rho = the value of parameter rho
%% options.index_alpha = the index set alpha 
%% options.index_beta = the index set beta 
%% options.r_indexJ = the cardinality of the index set indexJ
%% options.AJ = the matrix AJ 
%% options.same_indexJ = 1, indexJ has not changed since the last
%%                          SSN iteration
%%                     = 0, otherwise
%% options.num_T = 1, compute the matrix T
%%               = 0, does not compute it
%% options.gamISzero = 1, if const_gam = 0
%%                   = 0, otherwise
%% options.Palpha = submatrix of P by keeping all the rows in alpha
%% options.Pbeta = submatrix of P by keeping all the rows in beta
%% options.Pgamma = submatrix of P by keeping all the rows in gamma
%% options.AJTAJ = AJ'*AJ
%% Output:
%% invDrhs = the vector D^{-1}*rhs
%% T = the matrix T
%% options.AJTAJ = AJ'*AJ
%% Copyright (c) 2022 by Can Wu, Ying Cui, Donghui Li, Defeng Sun
%% For more details, please see the paper: 
%% Convex and Nonconvex Risk-based Linear Regression at Scale.
%%=========================================================================
function [invDrhs,T,options] = invDvec_MM(rhs,flag_case,options,computeT)
sigma = options.sigma;
eta = options.eta;
rho = options.rho;

index_alpha = options.index_alpha;
index_beta = options.index_beta;
r_indexJ = options.r_indexJ;
AJ = options.AJ;
same_indexJ = options.same_indexJ;
num_T = options.num_T;
gamISzero = options.gamISzero;

Palpha = options.Palpha;
Pbeta = options.Pbeta;
Pgamma = options.Pgamma;
T = [];

delta = 1/sigma + 1/eta;
if (computeT == 1)
    if (same_indexJ == 1) && (num_T > 0) && (size(options.AJTAJ,1) == r_indexJ)
        AJTAJ = options.AJTAJ;
    else
        AJTAJ = AJ'*AJ;
    end
    options.AJTAJ = AJTAJ;
end
t_tmp =  eta/(eta+sigma);
num2 = eta/sigma;

switch(flag_case)
    case 1
        invDrhs = eta*rhs;
        if (computeT == 1) && (r_indexJ > 0)
            T = rho*speye(r_indexJ) + eta*AJTAJ;
        end
        return;
    case {2.1, 2.2}
        if (flag_case == 2.1)
            Pbet_gam = [Pbeta;Pgamma];
            invDrhs = (1/delta)*(rhs + num2*((Pbet_gam*rhs)'*Pbet_gam)');
            if (computeT == 1) && (r_indexJ > 0)
                Pbet_gamAJ = Pbet_gam*AJ;
                T = rho*speye(r_indexJ) + (1/delta)*(AJTAJ + num2*(Pbet_gamAJ'*Pbet_gamAJ));
            end
            return;
        else
            invDrhs = eta*(rhs - t_tmp*((Palpha*rhs)'*Palpha)');
            if (computeT == 1) && (r_indexJ > 0)
                C_tmp = Palpha*AJ;
                T = rho*speye(r_indexJ) + eta*(AJTAJ-t_tmp*(C_tmp'*C_tmp));
            end
            return;
        end
    case {3.1, 3.2}
        const3 = index_alpha+index_beta;
        Palp_bet = [Palpha;Pbeta];
        if (flag_case == 3.1)
            V31 = [ones(1,const3)*Palp_bet; Palp_bet];
            rhs32 = Palp_bet*rhs;
            rhs31 = ones(1,const3)*rhs32;
            num_d = t_tmp/const3;
            d31 = num_d*rhs31;
            d3 =[d31;num2*(d31*ones(const3,1)-rhs32)];
            invDrhs = (1/delta)*(rhs-(d3'*V31)');
            if (computeT == 1) && (r_indexJ > 0)
                V31AJ = V31*AJ;
                R3 = Palp_bet*AJ;
                d_tmp = num_d*(ones(1,const3)*R3);
                Tinv_tmp = [d_tmp;num2*(ones(const3,1)*d_tmp-R3)];
                T = rho*speye(r_indexJ) + (1/delta)*( AJTAJ-(V31AJ)'*Tinv_tmp );
            end
            return;
        else
            ePalp_bet = ones(1,const3)*Palp_bet;
            if gamISzero
                invDrhs = eta*(rhs - t_tmp* (1/const3)*(ePalp_bet*rhs)*ePalp_bet');
            else
                invDrhs = eta*(rhs - t_tmp*( (1/const3)*(ePalp_bet*rhs)*ePalp_bet'+ ((Pgamma*rhs)'*Pgamma)'));
            end
            if (computeT == 1) && (r_indexJ > 0)
                if gamISzero
                    C_tmp = (1/sqrt(const3))*(ePalp_bet*AJ);
                else
                    C_tmp = [(1/sqrt(const3))*(ePalp_bet*AJ); Pgamma*AJ ];
                end
                T = rho*speye(r_indexJ)+eta*(AJTAJ-t_tmp*(C_tmp'*C_tmp));
            end
            return;
        end
    case {4.1, 4.2}
        if (flag_case == 4.1)
            V41 = [Pbeta;ones(1,index_beta)*Pbeta];
            rhs41 = Pbeta*rhs;
            rhs42 = ones(1,index_beta)*rhs41;
            num_d = (t_tmp/index_beta);
            d42 = num_d*rhs42;
            d4 = [num2*(d42*ones(index_beta,1)-rhs41);d42];
            invDrhs = (1/delta)*(rhs-(d4'*V41)');
            if (computeT == 1) && (r_indexJ > 0)
                R4 = Pbeta*AJ;
                d_tmp = num_d*(ones(1,index_beta)*R4);
                Tinv_tmp = [num2*(ones(index_beta,1)*d_tmp-R4);d_tmp];
                T = rho*speye(r_indexJ) + (1/delta)*( AJTAJ-((V41*AJ)'*Tinv_tmp) );
            end
            return;
        else
            ePbeta = ones(1,index_beta)*Pbeta;
            if gamISzero
                invDrhs = eta*( rhs-t_tmp*( ((Palpha*rhs)'*Palpha)'+((ePbeta*rhs)/index_beta)*ePbeta'));
            else
                invDrhs = eta*( rhs-t_tmp*( ((Palpha*rhs)'*Palpha)'+((ePbeta*rhs)/index_beta)*ePbeta'+((Pgamma*rhs)'*Pgamma)'));
            end
            if (computeT == 1) && (r_indexJ > 0)
                if gamISzero
                    C_tmp = [Palpha*AJ; (1/sqrt(index_beta))*(ePbeta*AJ)];
                else
                    C_tmp = [Palpha*AJ; (1/sqrt(index_beta))*(ePbeta)*AJ; Pgamma*AJ];
                end
                T = rho*speye(r_indexJ) + eta*(AJTAJ-t_tmp*(C_tmp'*C_tmp));
            end
            return;
        end
    case 5
        invDrhs = (1/delta)*rhs;
        if (computeT == 1) && (r_indexJ > 0)
            T = rho*speye(r_indexJ)+(1/delta)*AJTAJ;
        end
        return;
end


end