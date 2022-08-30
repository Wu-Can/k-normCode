function [invDrhs,T,options] = invDvec(rhs,flag_case,options,computeT)
%% This function is to calculate the vector D^{-1}*rhs
%% and to calculate the matrix T = sigma^{-1}*I_s + AJ^{T}*D^{-1}*AJ
%% if computeT = 1.

sigma = options.sigma;
tau = options.tau;
index_alpha = options.index_alpha;
index_beta = options.index_beta;
index_gamma = options.index_gamma;
r_indexJ = options.r_indexJ;
AJ = options.AJ;
same_indexJ = options.same_indexJ;
num_T = options.num_T;

Palpha = options.Palpha;
Pbeta = options.Pbeta;
Pgamma = options.Pgamma;
T = [];

delta = sigma + tau/sigma;
if (computeT == 1)
    if (same_indexJ == 1) && (num_T > 0) && (size(options.AJTAJ,1) == r_indexJ)
        AJTAJ = options.AJTAJ;
    else
        AJTAJ = AJ'*AJ;
    end
    options.AJTAJ = AJTAJ;
end

theta2 = 1;
t_tmp =  (sigma^2)/(tau+sigma^2);
num2 = sigma^2/tau;

switch(flag_case)
    case 1
        invDrhs = (sigma/tau)*rhs;
        if (computeT == 1) && (r_indexJ > 0)
            T = (1/sigma)*speye(r_indexJ)+(sigma/tau)*AJTAJ;
        end
        return;
    case 2
        if index_alpha >= theta2*index_beta
            invDrhs = (1/delta)*(rhs + ((sigma^2)/tau)*((Pbeta*rhs)'*Pbeta)');
            if (computeT == 1) && (r_indexJ > 0)
                PbetaAJ = Pbeta*AJ;
                T = (1/sigma)*speye(r_indexJ) + (1/delta)*(AJTAJ + ((sigma^2)/tau)*(PbetaAJ'*PbetaAJ));
            end
            options.flag_small2 = 1;
            return;
        else
            invDrhs = (sigma/tau)*(rhs - t_tmp*((Palpha*rhs)'*Palpha)');
            if (computeT == 1) && (r_indexJ > 0)
                C_tmp = Palpha*AJ;
                T = (1/sigma)*speye(r_indexJ)+(sigma/tau)*(AJTAJ-t_tmp*(C_tmp'*C_tmp));
            end
            options.flag_small2 = 0;
            return;
        end
    case {3.1, 3.2}
        if index_alpha <= theta2*index_beta
            V31 = [ones(1,index_alpha)*Palpha; Palpha];
            rhs32 = Palpha*rhs;
            rhs31 = ones(1,index_alpha)*rhs32;
            num_d = (t_tmp/index_alpha);
            d31 = num_d*rhs31;
            d3 =[d31;num2*(d31*ones(index_alpha,1)-rhs32)];
            invDrhs = (1/delta)*(rhs-(d3'*V31)');
            if (computeT == 1) && (r_indexJ > 0)
                V31AJ = V31*AJ;
                R3 = Palpha*AJ;
                d_tmp = num_d*(ones(1,index_alpha)*R3);
                Tinv_tmp = [d_tmp;num2*(ones(index_alpha,1)*d_tmp-R3)];
                T = (1/sigma)*speye(r_indexJ) + (1/delta)*( AJTAJ-(V31AJ)'*Tinv_tmp );
            end
            options.flag_small2 = 1;
            return;
        else
            ePalpha = ones(1,index_alpha)*Palpha;
            if flag_case == 3.1
                invDrhs = (sigma/tau)*(rhs - t_tmp* (1/index_alpha)*(ePalpha*rhs)*ePalpha');
            else
                invDrhs = (sigma/tau)*(rhs - t_tmp*( (1/index_alpha)*(ePalpha*rhs)*ePalpha'+ ((Pbeta*rhs)'*Pbeta)'));
            end
            if (computeT == 1) && (r_indexJ > 0)
                if flag_case == 3.1
                    C_tmp = (1/sqrt(index_alpha))*ePalpha*AJ;
                else
                    C_tmp = [(1/sqrt(index_alpha))*ePalpha*AJ; Pbeta*AJ ];
                end
                T = (1/sigma)*speye(r_indexJ)+(sigma/tau)*(AJTAJ-t_tmp*(C_tmp'*C_tmp));
            end
            options.flag_small2 = 0;
            return;
        end
    case {4.1, 4.2}
        if index_beta < theta2*(index_alpha+index_gamma)
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
                T = (1/sigma)*speye(r_indexJ) + (1/delta)*( AJTAJ-((V41*AJ)'*Tinv_tmp) );
            end
            options.flag_small2 = 1;
            return;
        else
            ePbeta = ones(1,index_beta)*Pbeta;
            if (flag_case == 4.1)
                invDrhs = (sigma/tau)*( rhs-t_tmp*( ((Palpha*rhs)'*Palpha)'+((ePbeta*rhs)/index_beta)*ePbeta'));
            else
                invDrhs = (sigma/tau)*( rhs-t_tmp*( ((Palpha*rhs)'*Palpha)'+((ePbeta*rhs)/index_beta)*ePbeta'+((Pgamma*rhs)'*Pgamma)'));
            end
            if (computeT == 1) && (r_indexJ > 0)
                if (flag_case == 4.1)
                    C_tmp = [Palpha*AJ; (1/sqrt(index_beta))*(ePbeta)*AJ];
                else
                    C_tmp = [Palpha*AJ; (1/sqrt(index_beta))*(ePbeta)*AJ; Pgamma*AJ];
                end
                T = (1/sigma)*speye(r_indexJ)+(sigma/tau)*(AJTAJ-t_tmp*(C_tmp'*C_tmp));
            end
            options.flag_small2 = 0;
            return;
        end
    case 5
        invDrhs = (1/delta)*rhs;
        if (computeT == 1) && (r_indexJ > 0)
            T = (1/sigma)*speye(r_indexJ)+(1/delta)*AJTAJ;
        end
        return;
end

end